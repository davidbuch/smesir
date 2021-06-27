
smesir_forecast <- function(J,Jf,xi_samps,sigma2_samps,prior = NULL){
  # 3. Check that the prior is a list of the valid form  
  if(is.null(prior)){
    prior <- list(ell = J/5, V0 = c(10,10,1), IGSR = matrix(rep(2,6), nrow = 3, ncol = 2))
  }else{
    if(K == 1){
      if(!is.numeric(prior[["ell"]]) || length(prior[["ell"]]) != 1 || prior[["ell"]] <= 0){
        stop("prior[['ell']] must be a positive scalar")
      }
      if(!is.numeric(prior[["V0"]]) || length(prior[["V0"]]) != 2 || any(prior[["ell"]] <= 0)){
        stop("prior[['V0']] must be a length 2 positive numeric vector (single-region style)")
      }
      if(!is.numeric(prior[["IGSR"]]) || length(prior[["IGSR"]]) != 2 || any(prior[["IGSR"]] <= 0)){
        stop("prior[['IGSR']] must be a length 2 positive numeric vector (single-region style)")
      }
      prior[["V0"]] <- c(prior[["V0"]],1) # last value placeholder
      prior[["IGSR"]] <- matrix(rep(prior[["IGSR"]],3), nrow = 3, ncol = 2,  byrow = TRUE)
    }else{
      if(!is.numeric(prior[["ell"]]) || length(prior[["ell"]]) != 1 || prior[["ell"]] <= 0){
        stop("prior[['ell']] must be a positive scalar")
      }
      if(!is.numeric(prior[["V0"]]) || length(prior[["V0"]]) != 3 || any(prior[["ell"]] <= 0)){
        stop("prior[['V0']] must be a length 3 positive numeric vector")
      }
      if(!is.numeric(prior[["IGSR"]]) || dim(prior[["IGSR"]]) != c(3,2) || any(prior[["IGSR"]] <= 0)){
        stop("prior[['IGSR']] must be a 3x2 positive numeric matrix")
      }
    }
  }
  ell <- prior[["ell"]]
  V0 <- prior[["V0"]]
  IGSR <- prior[["IGSR"]]
  
  S_fit <- exp(-as.matrix(dist(1:J, diag = TRUE, upper = TRUE)/ell)^2)
  S_fit <- scale(S_fit, center = TRUE, scale = FALSE) # Project out the intercept
  kernel_decomp <- eigen(S_fit,symmetric = TRUE)
  rfit <- match(TRUE, cumsum(kernel_decomp$values)/sum(kernel_decomp$values) > .99)
  ebasis <- kernel_decomp$vectors[,1:rfit]
  lambda <- kernel_decomp$values[1:rfit]

  
  Jt = J + Jf;
  
  S <- exp(-as.matrix(dist(1:Jt, diag = TRUE, upper = TRUE)/ell)^2)
  S <- Re(scale(S, center = TRUE, scale = FALSE)) # Project out the intercept

  
  Li <- diag(1/lambda)
  S_11 <- S[1:J,1:J]
  S_12 <- S[(J + 1):Jt,1:J]
  S_22 <- S[(J+1):Jt,(J+1):Jt]
  
  mu_transform <- S_12%*%ebasis%*%Li
  conditional_cov <-  S_22 - S_12%*%ebasis%*%Li%*%t(ebasis)%*%t(S_12)
  
  kernel_decomp <- eigen(conditional_cov,symmetric = TRUE)
  kernel_decomp$values <- sapply(kernel_decomp$values,function(x) max(x,0))
  rf <- match(TRUE, cumsum(kernel_decomp$values)/sum(kernel_decomp$values) > .99)
  ebasis_forecasting <- kernel_decomp$vectors[,1:rf]
  colnames(ebasis_forecasting) <- paste("GPf Basis Func.",1:rf)
  vscales_theta_forecasting <- kernel_decomp$values[1:rf]
  
  # apply mu_transform to samples of xi to get the predictive mean
  # draw independent gaussians with variances "vscales_theta_forecasting"
  # premultiply those coefficients by ebasis_forecasting and add them to 
  # the predictive means we have determined.
  # Thereby we will get a sample from the posterior
  
  forecast_expectations <- xi_samps%*%t(mu_transform)
  forecast_theta <- sapply(sigma2_samps, function(sigma2){rnorm(rf, sd = sqrt(sigma2*vscales_theta_forecasting))})
  forecast_noise <- t(ebasis_forecasting%*%forecast_theta)
  forecast_beta <- forecast_expectations + forecast_noise
  
  return(forecast_beta)
}