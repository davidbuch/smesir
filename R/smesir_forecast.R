#' smesir_forecast: 
#' @description
#' Forecasts incidence events within the smeSIR model framework. A companion
#' function to \code{smesir}.
#' 
#' @param Jf Integer length of the forecast time horizon
#' @param model_fit A list returned from the \code{smesir} function
#' @param new_x An list containing covariate data with which to forecast, as needed by
#' by the model fit in the antecedent "model_fit" call.
#' @export
smesir_forecast <- function(Jf, model_fit, new_x = NULL, new_vaccinations = NULL){
  if(Jf < 1 || (Jf %% 1 != 0)) stop("Jf (forecast time horizon) must be an integer greater than 0.")
  
  ## Unpack the model_fit object and preliminary tests
  prior <- model_fit[["prior"]]
  formula <- model_fit[["formula"]]
  xi_samps <- model_fit[["samples"]][["Xi"]]
  
  N <- model_fit$epi_params$region_populations
  T_1 <- model_fit$epi_params$outbreak_times
  gamma <- 1/model_fit$epi_params$mean_removal_time
  psi <- model_fit$epi_params$incidence_probabilities
    
  design_matrices <- model_fit[["design_matrices"]]
  
  ell <- prior[["ell"]]
  V0 <- prior[["V0"]]
  IGSR <- prior[["IGSR"]]
  
  K <- length(design_matrices)
  Jo <- nrow(design_matrices[[1]])
  nsamps <- nrow(xi_samps)
  
  # if the user would like to specify, extend this to include "new vaccinations"
  # otherwise maybe extrapolate trend from last 4 weeks
  # Also need to switch event forecasts to NB with dispersion samples.
  if(is.null(new_vaccinations)){
    vaccinations <- rbind(model_fit$vaccinations,matrix(0, nrow = Jf, ncol = K))
  }else if(any(dim(new_vaccinations) != c(Jf,K))){
    stop("Dimensions of 'new_vaccinations' must equal 'c(Jf,K)'")
  }else{
    vaccinations <- rbind(model_fit$vaccinations,new_vaccinations)
  }
  
  covariate_names <- all.vars(formula[[3]]) # this is really nice because it throws errors for one sided formulas
  if(length(covariate_names) > 0){
    for(cname in covariate_names){
      if(!(cname %in% names(new_x))) stop(paste0("Required covariate '",cname,"' not provided"))
      if(K == 1){
        if(length(new_x[[cname]]) != Jf) stop("Rows in 'new_x' must equal 'Jf'.")
      }else{
        if(is.null(dim(new_x[[cname]])) || any(dim(new_x[[cname]]) != c(Jf,K))) stop(paste0("Covariate '", cname,"' in new_x has incorrect dimensions"))
      }
    }
  }
  response_name <- all.vars(formula[[2]])
  if(length(response_name) != 1) stop("Error in argument 'formula': Response must be univariate.")
  if(!attr(terms(formula),"intercept")) stop("Error in argument 'formula': smeSIR model must have an intercept.")
  
  sr_style <- model_fit[["sr_style"]]
  if(!sr_style){
    sigma2_samps <- model_fit[["samples"]][["V"]][,3]
  }else{
    sigma2_samps <- model_fit[["samples"]][["V"]]
  }
  
  ## Create forecasting design matrices
  # Recreate GP covariance matrix from the model fit
  S_fit <- exp(-as.matrix(dist(1:Jo, diag = TRUE, upper = TRUE)/ell)^2)
  S_fit <- scale(S_fit, center = TRUE, scale = FALSE) # Project out the intercept
  kernel_decomp <- eigen(S_fit,symmetric = TRUE)
  rfit <- match(TRUE, cumsum(kernel_decomp$values)/sum(kernel_decomp$values) > .99)
  ebasis <- kernel_decomp$vectors[,1:rfit]
  lambda <- kernel_decomp$values[1:rfit]
  
  # Create GP covariance matrix for the complete model
  Jt = Jo + Jf
  S <- exp(-as.matrix(dist(1:Jt, diag = TRUE, upper = TRUE)/ell)^2)
  S <- Re(scale(S, center = TRUE, scale = FALSE)) # Project out the intercept
  
  Li <- diag(1/lambda)
  S_11 <- S[1:Jo,1:Jo]
  S_12 <- S[(Jo + 1):Jt,1:Jo]
  S_22 <- S[(Jo+1):Jt,(Jo+1):Jt]
  
  mu_transform <- S_12%*%ebasis%*%Li
  conditional_cov <-  S_22 - S_12%*%ebasis%*%Li%*%t(ebasis)%*%t(S_12)
  
  kernel_decomp <- eigen(conditional_cov,symmetric = TRUE)
  kernel_decomp$values <- sapply(kernel_decomp$values,function(x) max(x,0))
  rf <- match(TRUE, cumsum(kernel_decomp$values)/sum(kernel_decomp$values) > .99)
  ebasis_forecasting <- kernel_decomp$vectors[,1:rf]
  colnames(ebasis_forecasting) <- paste("GP Basis Func.",1:rf)
  vscales_theta_forecasting <- kernel_decomp$values[1:rf]
  
  design_matrices_forecast <- list()
  if(K == 1){
    if(is.null(dim(new_x))){
      stop("When only forecasting for one region, 'new_x' must be a data frame.")
    }
    covariate_dframe <- list()
    for(cname in covariate_names){
      covariate_dframe[[cname]] <- new_x[,cname]
    }
    covariate_dframe <- as.data.frame(covariate_dframe)
    design_mat <- model.matrix(as.formula(paste("~", as.character(formula)[3])), covariate_dframe)
    if(length(covariate_names) == 0){
      design_mat <- matrix(1, nrow = Jf, ncol = 1)
      colnames(design_mat) <- "(Intercept)"
    }
    design_matrices_forecast[[1]] <- cbind(design_mat,ebasis_forecasting)
  }else{
    for(k in 1:K){
      covariate_dframe <- list()
      for(cname in covariate_names){
        covariate_dframe[[cname]] <- new_x[[cname]][,k]
      }
      covariate_dframe <- as.data.frame(covariate_dframe)
      design_mat <- model.matrix(as.formula(paste( "~", as.character(formula)[3])), covariate_dframe)
      if(length(covariate_names) == 0){
        design_mat <- matrix(1, nrow = Jf, ncol = 1)
        colnames(design_mat) <- "(Intercept)"
      }
      design_matrices_forecast[[k]] <- cbind(design_mat,ebasis_forecasting) ## attach eigenbasis
    }
  }
  P <- ncol(design_matrices[[1]]) # number of "predictors" (including intercept and GP bases)
  nterms <- P - rfit - 1 # not the same as the number of covariates (e.g., factor expansions)
  predictor_names <- colnames(design_matrices[[1]])
  vparam_names <- c("Variance(Intercept)", "Variance(Covariate Coeffs.)", "Variance(GP Random Effect)")
  
  
  ## Transform coefficient samples into forecast samples
  if(K == 1){
    forecast_beta <- array(dim = c(Jf,nsamps))
    beta_samps <- array(dim = c(Jo + Jf, nsamps))    
    
    forecast_expectations <- mu_transform %*% t(xi_samps[,(1+nterms+1):P])
    forecast_theta <- sapply(sigma2_samps, function(sigma2){rnorm(rf, sd = sqrt(sigma2*vscales_theta_forecasting))})
    forecast_noise <- ebasis_forecasting%*%forecast_theta
    forecast_beta <- forecast_expectations + forecast_noise + as.matrix(design_matrices_forecast[[1]][,1:(1+nterms)])%*%t(xi_samps[,1:(1+nterms)])
    beta_samps <- rbind(design_matrices[[1]] %*% t(xi_samps),forecast_beta)
    beta_samps <- pmax(beta_samps,0)
    
    event_samps <- array(dim = c(Jo + Jf, nsamps))
    event_CI <- array(dim = c(Jo + Jf, 2))
    
    event_comp <- function(beta_samp,iip_samp){
      rpois(length(beta_samp), solve_events(solve_infections(beta_samp,
                                                            gamma, T_1,
                                                            iip_samp, N, 
                                                            vaccinations),psi))
    }
    for(s in 1:nsamps){
      event_samps[,s] <- event_comp(beta_samps[,s],model_fit$samples$IIP[s])
    }
    event_CI <- t(apply(event_samps,1,function(x) quantile(x, c(0.025,0.975))))
  }else{
    forecast_beta <- array(dim = c(Jf,nsamps,K))
    beta_samps <- array(dim = c(Jo + Jf, nsamps,K))
    for(k in 1:K){
      forecast_expectations <- mu_transform %*% t(xi_samps[,(1+nterms+1):P,k])
      forecast_theta <- sapply(sigma2_samps, function(sigma2){rnorm(rf, sd = sqrt(sigma2*vscales_theta_forecasting))})
      forecast_noise <- ebasis_forecasting%*%forecast_theta
      forecast_beta[,,k] <- forecast_expectations + forecast_noise + as.matrix(design_matrices_forecast[[k]][,1:(1+nterms)])%*%t(xi_samps[,1:(1+nterms),k])
      beta_samps[,,k] <- rbind(design_matrices[[k]] %*% t(xi_samps[,,k]),forecast_beta[,,k])
      beta_samps[,,k] <- pmax(beta_samps[,,k],0)
      # plot(rowMeans(beta_samps[,,k]), type = "l", ylim = c(0,3), main = k)
      # plot_confint(t(apply(beta_samps[,,k],1,function(x) quantile(x,c(0.025,0.975)))),density=15,col = "blue")
    }
    
    event_samps <- array(dim = c(Jo + Jf, nsamps,K))
    event_CI <- array(dim = c(Jo + Jf, 2,K))
    for(k in 1:K){
      event_comp <- function(beta_samp,iip_samp){
        rpois(length(beta_samp), solve_events(solve_infections(beta_samp,
                                                              gamma, T_1[k],
                                                              iip_samp, N[k],
                                                              vaccinations[,k]),psi))
      }
      for(s in 1:nsamps){
        event_samps[,s,k] <- event_comp(beta_samps[,s,k],model_fit$samples$IIP[s,k])
      }
      event_CI[,,k] <- t(apply(event_samps[,,k],1,function(x) quantile(x, c(0.025,0.975))))
    }
  }
  forecast_dat <- list(samples = event_samps, confints = event_CI)
  
  return(forecast_dat)
}