#' smesir: 
#' @description
#' Forecasts incidence events 
#' 
#' @param Jf Integer length of the forecast time horizon
#' @param model_fit A list returned from the \code{smesir} function
#' @param new_x An optional list containing covariate data to be used in forecasting
smesir_forecast <- function(Jf, model_fit, new_x = NULL, region_populations, outbreak_times, 
                            mean_removal_time, incidence_probabilities, dispersion,
                            initial_impulses = 1/region_populations){
  # unpack the model_fit object
  prior <- model_fit[["prior"]]
  formula <- model_fit[["formula"]]
  xi_samps <- model_fit[["samples"]][["Xi"]]
  
  ell <- prior[["ell"]]
  V0 <- prior[["V0"]]
  IGSR <- prior[["IGSR"]]
  
  covariate_names <- all.vars(formula[[3]]) # this is really nice because it throws errors for one sided formulas
  response_name <- all.vars(formula[[2]])
  if(length(response_name) != 1) stop("Error in argument 'formula': Response must be univariate.")
  if(!attr(terms(formula),"intercept")) stop("Error in argument 'formula': smeSIR model must have an intercept.")
  
  
  if(length(dim(xi_samps)) == 3){
    K <- dim(xi_samps)[3]
  }else{
    K <- 1
  }
  sr_style <- model_fit[["sr_style"]]
  if(!sr_style){
    sigma2_samps <- model_fit[["samples"]][["V"]][,3]
  }else{
    sigma2_samps <- model_fit[["samples"]][["V"]]
  }
  
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
  
  design_matrices_forecast <- list()
  if(K == 1){
    covariate_dframe <- list()
    for(cname in covariate_names){
      covariate_dframe[[cname]] <- new_x[,cname]
    }
    covariate_dframe <- as.data.frame(covariate_dframe)
    design_mat <- model.matrix(as.formula(paste("~", as.character(formula)[3])), covariate_dframe)
    if(length(covariate_names) == 0){
      design_mat <- matrix(1, nrow = J, ncol = 1)
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
        design_mat <- matrix(1, nrow = J, ncol = 1)
        colnames(design_mat) <- "(Intercept)"
      }
      design_matrices_forecast[[k]] <- cbind(design_mat,ebasis_forecasting) ## attach eigenbasis
    }
  }
  ncov <- ncol(design_mat) - 1
  P <- ncol(design_matrices_forecast[[1]]) # number of "predictors"
  predictor_names <- colnames(design_matrices_forecast[[1]])
  vparam_names <- c("Variance(Intercept)", "Variance(Covariate Coeffs.)", "Variance(GP Random Effect)")
  # names(design_matrices_forecast) <- region_names
   

  # apply mu_transform to samples of xi to get the predictive mean
  # draw independent gaussians with variances "vscales_theta_forecasting"
  # premultiply those coefficients by ebasis_forecasting and add them to 
  # the predictive means we have determined.
  # Thereby we will get a sample from the posterior
  forecast_beta <- array(dim = c(nsamps,Jf,K))
  for(k in 1:K){
    forecast_expectations <- xi_samps[,(1+ncov+1):P,k]%*%t(mu_transform)
    forecast_theta <- sapply(sigma2_samps, function(sigma2){rnorm(rf, sd = sqrt(sigma2*vscales_theta_forecasting))})
    forecast_noise <- t(ebasis_forecasting%*%forecast_theta)
    #print(dim(forecast_expectations))
    #print(dim(forecast_noise))
    forecast_beta[,,k] <- forecast_expectations + forecast_noise + xi_samps[,1:(1+ncov),k]%*%t(as.matrix(design_matrices_forecast[[k]][,1:(1+ncov)]))
  }

  
  deaths_comp <- function(beta_samp){
    rpois(length(beta_samp),solve_events(solve_infections(beta_samp,
                                                          gamma, outbreak_times[k],
                                                          1/region_pops[k], region_pops[k]),psi))
  }
  
  beta_samps <- array(dim = c(nsamps,J,K))
  for(k in 1:K){
    beta_samps[,,k] <- mfit$samples$Xi[,,k] %*% t(mfit$design_matrices[[k]])
    matplot(t(apply(cbind(beta_samps[,,k],forecast_beta[,,k]),2,function(x) quantile(x,c(0.025,0.5,0.975)))), type = "l")
  }
  
  death_samps <- array(dim = c(J + Jf,nsamps,K))
  deaths_CI <- array(dim = c(J + Jf,2,K))
  for(k in 1:K){
    deaths_comp <- function(beta_samp){
      if(dispersion > 0){
      rnbinom(length(beta_samp),size = 1/dispersion, mu = solve_events(solve_infections(beta_samp,
                                                            1/mean_removal_time, outbreak_times[k],
                                                            initial_impulses[k], region_populations[k]),
                                           incidence_probabilities))
      }else if(dispersion == 0){
        rpois(length(beta_samp),solve_events(solve_infections(beta_samp,
                                                              1/mean_removal_time, outbreak_times[k],
                                                              initial_impulses[k], region_populations[k]),
                                             incidence_probabilities))
      }else{stop("Invalid dispersion parameter")}
    }
    #print(dim(apply(cbind(beta_samps[,,k],forecast_beta[,,k]),1,deaths_comp)))
    death_samps[,,k] <- apply(cbind(beta_samps[,,k],forecast_beta[,,k]),1,deaths_comp)
    #plot(rowMeans(death_samps[,,k]),type = "l")
    #print(dim(apply(death_samps[,,k],1,function(x) quantile(x,c(0.025,0.975)))))
    deaths_CI[,,k] <- t(apply(death_samps[,,k],1,function(x) quantile(x,c(0.025,0.975))))
  }
  forecast_dat <- list(samples = death_samps, summary = deaths_CI)
  
  return(forecast_dat)
}