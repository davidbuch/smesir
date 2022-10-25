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
smesir_forecast <- function(Jf, model_fit, new_x = NULL, new_vaccinations = NULL, fixed_dispersion = NULL){
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
  
  r <- model_fit$n_basis
  
  K <- length(design_matrices)
  Jo <- nrow(design_matrices[[1]])
  P <- ncol(design_matrices[[1]])
  Jt <- Jo + Jf
  nsamps <- nrow(xi_samps)
  
  # if the user would like to specify, extend this to include "new vaccinations"
  # otherwise assume no more vaccinations
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
  
  forecast_covariates <- list()
  if(K == 1){
    if(!is.null(new_x) && is.null(dim(new_x))){
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
    forecast_covariates[[1]] <- design_mat
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
      forecast_covariates[[k]] <- design_mat
    }
  }
 # number of "predictors" (including intercept and GP bases)
  nterms <- P - r - 1 # not the same as the number of covariates (e.g., factor expansions)
  predictor_names <- colnames(design_matrices[[1]])
  vparam_names <- c("Variance(Intercept)", "Variance(Covariate Coeffs.)", "Variance(GP Random Effect)")
  
  
  S <- exp(-as.matrix(dist(1:Jt, diag = TRUE, upper = TRUE)/ell)^2)

  S_11 <- S[1:Jo,1:Jo] # After projection, should be the same as "smesir" created
  S_21 <- S[(Jo + 1):Jt,1:Jo]
  S_22 <- S[(Jo+1):Jt,(Jo+1):Jt]
  
  eS_11 <- eigen(S_11)
  ebasis <- eS_11$vectors[,1:r]
  Li <- diag(1/eS_11$values[1:r])
  
  conditional_mean_transform <- S_21%*%ebasis%*%Li
  eS2g1 <- eigen(S_22 - S_21%*%ebasis%*%Li%*%t(ebasis)%*%t(S_21))
  rf <- min(Jf,r)
  conditional_cov_transform <- eS2g1$vectors[,1:rf] %*% diag(sqrt(eS2g1$values[1:rf]))
  
  
  ## Transform coefficient samples into forecast samples
  if(K == 1){
    forecast_expectations <- conditional_mean_transform %*% t(xi_samps[,(P-r+1):P])
    forecast_gp <- forecast_expectations + conditional_cov_transform %*% matrix(rnorm(rf*nsamps), nrow = rf, ncol = nsamps)
    forecast_beta <- forecast_gp + as.matrix(forecast_covariates[[1]])%*%t(xi_samps[,1:(1+nterms)])
    beta_samps <- exp(rbind(design_matrices[[1]] %*% t(xi_samps),forecast_beta))
    #beta_samps <- pmax(beta_samps,0)
#    plot(t(apply(beta_samps/gamma,1,function(x) quantile(x,0.5)), type = "l", ylim = c(0,max(apply(beta_samps/gamma,1,function(x) quantile(x,c(0.025,0.975))))), main = k)
#    plot_confint(t(apply(beta_samps/gamma,1,function(x) quantile(x,c(0.025,0.975)))),density=15,col = "blue")
    
    expected_samps <- array(dim = c(Jo + Jf, nsamps))
    event_samps <- array(dim = c(Jo + Jf, nsamps))
    vulnerable_samps <- array(dim = c(Jo + Jf, nsamps))
    event_CI <- array(dim = c(Jo + Jf, 3))
    expected_CI <- array(dim = c(Jo + Jf, 3))
    vulnerable_CI <- array(dim = c(Jo + Jf, 3))
    
    expected_comp <- function(beta_samp,iip_samp,disp_samp){
      disp_samp <- ifelse(is.null(fixed_dispersion), disp_samp, fixed_dispersion)
      solve_events(solve_infections(beta_samp,
                                    gamma, T_1,
                                    iip_samp, N,
                                    vaccinations),psi)
    }
    event_comp <- function(beta_samp,iip_samp,disp_samp){
      disp_samp <- ifelse(is.null(fixed_dispersion), disp_samp, fixed_dispersion)
      rnbinom(length(beta_samp),
              size = 1/disp_samp,
              mu = solve_events(solve_infections(beta_samp,gamma, T_1, iip_samp, N, vaccinations),psi))
    }
    vulnerable_comp <- function(beta_samp,iip_samp){
      solve_susceptible(beta_samp, gamma, T_1,iip_samp, N,vaccinations)
    }
    for(s in 1:nsamps){
      expected_samps[,s] <- expected_comp(beta_samps[,s],model_fit$samples$IIP[s], model_fit$samples$DISP[s])
      event_samps[,s] <- event_comp(beta_samps[,s],model_fit$samples$IIP[s], model_fit$samples$DISP[s])
      vulnerable_samps[,s] <- vulnerable_comp(beta_samps[,s],model_fit$samples$IIP[s])
    }
    expected_CI <- t(apply(expected_samps,1,function(x) quantile(x, c(0.025,0.5,0.975))))
    event_CI <- t(apply(event_samps,1,function(x) quantile(x, c(0.025,0.5,0.975))))
    vulnerable_CI <- t(apply(vulnerable_samps,1,function(x) quantile(x, c(0.025,0.5,0.975))))
  }else{
    forecast_beta <- array(dim = c(Jf,nsamps,K))
    beta_samps <- array(dim = c(Jo + Jf, nsamps,K))
    for(k in 1:K){
      forecast_expectations <- conditional_mean_transform %*% t(xi_samps[,(P-r+1):P,k])
      forecast_gp <- forecast_expectations + conditional_cov_transform %*% matrix(rnorm(rf*nsamps), nrow = rf, ncol = nsamps)
      forecast_beta[,,k] <- forecast_gp + as.matrix(forecast_covariates[[k]])%*%t(xi_samps[,1:(1+nterms),k])
      beta_samps[,,k] <- exp(rbind(design_matrices[[k]] %*% t(xi_samps[,,k]),forecast_beta[,,k]))
      #beta_samps[,,k] <- pmax(beta_samps[,,k],0)
      
      #      plot(t(apply(beta_samps[,,k]/gamma,1,function(x) quantile(x,0.5))), type = "l", ylim = c(0,max(apply(beta_samps[,,k]/gamma,1,function(x) quantile(x,c(0.025,0.975))))), main = k)
      #      plot_confint(t(apply(beta_samps[,,k]/gamma,1,function(x) quantile(x,c(0.025,0.975)))),density=15,col = "blue")
    }
    
    expected_samps <- array(dim = c(Jo + Jf, nsamps,K))
    event_samps <- array(dim = c(Jo + Jf, nsamps,K))
    vulnerable_samps <- array(dim = c(Jo + Jf, nsamps,K))
    expected_CI <- array(dim = c(Jo + Jf, 3,K))
    event_CI <- array(dim = c(Jo + Jf, 3,K))
    vulnerable_CI <- array(dim = c(Jo + Jf, 3, K))
    for(k in 1:K){
      expected_comp <- function(beta_samp,iip_samp,disp_samp){
        disp_samp <- ifelse(is.null(fixed_dispersion), disp_samp, fixed_dispersion)
        solve_events(solve_infections(beta_samp, gamma, T_1[k],
                                      iip_samp, N[k],
                                      vaccinations[,k]),psi[,k])
      }
      event_comp <- function(beta_samp,iip_samp,disp_samp){
        disp_samp <- ifelse(is.null(fixed_dispersion), disp_samp, fixed_dispersion)
        rnbinom(length(beta_samp),
                size = 1/disp_samp,
                mu = solve_events(solve_infections(beta_samp, gamma, T_1[k],iip_samp,N[k],vaccinations[,k]),psi[,k]))
      }
      vulnerable_comp <- function(beta_samp,iip_samp){
        solve_susceptible(beta_samp, gamma, T_1[k],iip_samp, N[k],vaccinations[,k])
      }
      for(s in 1:nsamps){
        expected_samps[,s,k] <- expected_comp(beta_samps[,s,k],model_fit$samples$IIP[s,k],model_fit$samples$DISP[s,k])
        event_samps[,s,k] <- event_comp(beta_samps[,s,k],model_fit$samples$IIP[s,k],model_fit$samples$DISP[s,k])
        vulnerable_samps[,s,k] <- vulnerable_comp(beta_samps[,s,k],model_fit$samples$IIP[s,k])
      }
      expected_CI[,,k] <- t(apply(expected_samps[,,k],1,function(x) quantile(x, c(0.025,0.5,0.975))))
      event_CI[,,k] <- t(apply(event_samps[,,k],1,function(x) quantile(x, c(0.025,0.5,0.975))))
      vulnerable_CI[,,k] <- t(apply(vulnerable_samps[,,k],1,function(x) quantile(x, c(0.025,0.5,0.975))))
    }
  }
  forecast_dat <- list(expected_samples = expected_samps,
                       expected_confints = expected_CI,
                       event_samples = event_samps, 
                       event_confints = event_CI,
                       vulnerable_samples = vulnerable_samps,
                       vulnerable_confints = vulnerable_CI)
  
  return(forecast_dat)
}