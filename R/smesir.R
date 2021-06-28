#' smesir: Fitting Semiparametric Mixed-Effects SIR Models
#' 
#' @description
#' This function fits an SIR model with local, time-varying transmission 
#' rates to epidemic incidence data (deaths) observed in one or 
#' more regions. The transmission rate in each region is modeled
#' as a linear function of covariates plus a smooth temporal random 
#' effect drawn from a Gaussian Process distribution. The intercepts,
#' coefficients, and even the temporal random effects from the various
#' regions are assumed to be independent draws from global distributions
#' centered at "global" intercept, coefficient, and temporal random effect
#' values. 
#' 
#' @param formula an object of class "formula" 
#' @param data A length \code{P + 2} named list containing the data 
#' with which the model will be fit. Each name should match a name in
#' the provided \code{formula} argument
#' @param prior A length 2 named list containing \code{V0}, 3 variance hyperparameters for gaussian priors,
#' and \code{GSR}, 3 pairs of shape and rate hyperparameters for gamma priors
smesir <- function(formula, data, region_populations, outbreak_times, 
                   mean_removal_time, incidence_probabilities, 
                   initial_impulses = 1/region_populations, 
                   region_names = NULL, prior = NULL, inits = "random", 
                   chains = 4, iter = 2000, warmup = floor(iter / 2), thin = 1, 
                   min_adaptation_cycles = 10, min_samps_per_cycle = NULL, 
                   tempering_ratio = 0.2, quiet = TRUE, sr_style = NULL, 
                   seed = 12345){
  ## refer to the classic "lm" function docomentation for how this 
  ## function should be elaborated
  
  ## Sanity check all arguments
  
  # 0. check the existence and type of each variable
  # and region_pop and region_names should have length K
  # 1. check that each variable in "formula" is in the data list 
  # and also that formula includes an intercept.
  covariate_names <- all.vars(formula[[3]]) # this is really nice because it throws errors for one sided formulas
  response_name <- all.vars(formula[[2]])
  if(length(response_name) != 1) stop("Error in argument 'formula': Response must be univariate.")
  if(!attr(terms(formula),"intercept")) stop("Error in argument 'formula': smeSIR model must have an intercept.")
  # 2. check that elements of "data" list are the same size
  check2 <- FALSE
  JK <- dim(data[[1]])
  if(!is.null(dim(data[[1]]))){
    if(length(JK == 2)){ # no higer order tensors allowed
      check2 <- all(sapply(data,dim) == JK)
    }
    J <- JK[1]; K <- JK[2]
  }else{
    JK <- length(data[[1]])
    check2 <- all(sapply(data,length) == JK)
    J <- JK; K <- 1
  }
  if(is.null(sr_style)){ # user can override default sr_style behavior
    if(K == 1){
      sr_style = TRUE
    }else{
      sr_style = FALSE
    }
  }
  if(K == 1 && !sr_style) stop("sr_style set to FALSE but only one incidence sequence provided: Cannot fit multi-region style model to single stream of data.")
  if(!check2){
    stop("All elements of argument 'data' must be numeric matrices of equal size.")
  }
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
  
  ## We will need the eigendecomposition of the scaled GP covariance matrix
  ## when we construct our design matrix
  S_ell <- exp(-as.matrix(dist(1:J, diag = TRUE, upper = TRUE)/ell)^2)
  S_ell <- scale(S_ell, center = TRUE, scale = FALSE) # Project out the intercept
  kernel_decomp <- eigen(S_ell,symmetric = TRUE)
  r <- match(TRUE, cumsum(kernel_decomp$values)/sum(kernel_decomp$values) > .99)
  ebasis <- kernel_decomp$vectors[,1:r]
  colnames(ebasis) <- paste("GP Basis Func.",1:r)
  vscales_theta <- kernel_decomp$values[1:r]
  ## Construct Data Matrices for Model from user provided data list
  design_matrices <- list()
  response_matrix <- matrix(nrow = J, ncol = K)
  if(K == 1){
    covariate_dframe <- list()
    for(cname in covariate_names){
      covariate_dframe[[cname]] <- data[,cname]
    }
    covariate_dframe <- as.data.frame(covariate_dframe)
    design_mat <- model.matrix(as.formula(paste("~", as.character(formula)[3])), covariate_dframe)
    if(length(covariate_names) == 0){
      design_mat <- matrix(1, nrow = J, ncol = 1)
      colnames(design_mat) <- "(Intercept)"
    }
    design_matrices[[1]] <- cbind(design_mat,ebasis)
    #design_matrices[[1]] <- design_mat
    response_matrix[,1] <- data[,response_name]
  }else{
    for(k in 1:K){
      covariate_dframe <- list()
      for(cname in covariate_names){
        covariate_dframe[[cname]] <- data[[cname]][,k]
      }
      covariate_dframe <- as.data.frame(covariate_dframe)
      design_mat <- model.matrix(as.formula(paste( "~", as.character(formula)[3])), covariate_dframe)
      if(length(covariate_names) == 0){
        design_mat <- matrix(1, nrow = J, ncol = 1)
        colnames(design_mat) <- "(Intercept)"
      }
      design_matrices[[k]] <- cbind(design_mat,ebasis) ## attach eigenbasis
      #design_matrices[[k]] <- design_mat
      response_matrix[,k] <- data[[response_name]][,k]
    }
  }
  ncov <- ncol(design_mat) - 1
  P <- ncol(design_matrices[[1]]) # number of "predictors"
  predictor_names <- colnames(design_matrices[[1]])
  vparam_names <- c("Variance(Intercept)", "Variance(Covariate Coeffs.)", "Variance(GP Random Effect)")
  names(design_matrices) <- region_names
  colnames(response_matrix) <- region_names
  if(!sr_style){
    if(ncov == 0){
      V_used <- c(1,3)
    }else{
      V_used <- 1:3
    }
  }
  # when passing off to the MCMC functions, pass Y as a
  # matrix, whether or not K > 1. Similarly, pass XI as a matrix
  # and the Design_Matrices as a list (which is how I've already
  # constructed 'design_matrices' fortunately)
  set.seed(seed)
  if(!quiet){
    print("Reached the Sampling Function")
  }
  if(is.null(min_samps_per_cycle)){
    min_samps_per_cycle <- 10*P*P
  }
  MCMC_Output <- smesir_mcmc(response_matrix, design_matrices, vscales_theta, V0, IGSR, 1/mean_removal_time, outbreak_times, initial_impulses,
                             region_populations, incidence_probabilities, tempering_ratio, min_adaptation_cycles, min_samps_per_cycle, chains,iter,warmup,thin,sr_style,quiet) # last arg is sr_style flag
  
  ## do convergence diagnostics here
  nstore <- floor((iter - warmup)/thin)
  mcmc_diagnostics <- list()
  for(k in 1:K){
    # get them for Xi
    mcmc_diagnostics[[k]] <- matrix(nrow = P, ncol = 2)
    for(p in 1:P){
      samples_matrix <- matrix(nrow = nstore, ncol = chains)
      for(chn in 1:chains){
        source_row <- P*(k - 1) + p
        samples_matrix[,chn] <- MCMC_Output[[chn]][["Xi"]][source_row,]
      }
      mcmc_diagnostics[[k]][p,] <- convergence_diagnostics(samples_matrix)
    }
    rownames(mcmc_diagnostics[[k]]) <- predictor_names
    colnames(mcmc_diagnostics[[k]]) <- c("Rhat", "ESS")
  }
  if(!sr_style){
    # get them for Xi0 and all three variance params
    mcmc_diagnostics[[K + 1]] <- matrix(nrow = P+3, ncol = 2)
    for(p in 1:P){
      samples_matrix <- matrix(nrow = nstore, ncol = chains)
      for(chn in 1:chains){
        samples_matrix[,chn] <- MCMC_Output[[chn]][["Xi0"]][p,]
      }
      mcmc_diagnostics[[K + 1]][p,] <- convergence_diagnostics(samples_matrix)
    }
    for(vp in V_used){
      samples_matrix <- matrix(nrow = nstore, ncol = chains)
      for(chn in 1:chains){
        samples_matrix[,chn] <- MCMC_Output[[chn]][["V"]][vp,]
      }
      mcmc_diagnostics[[K + 1]][P + vp,] <- convergence_diagnostics(samples_matrix)
    }
    rownames(mcmc_diagnostics[[K + 1]]) <- c(predictor_names,vparam_names)
    colnames(mcmc_diagnostics[[K + 1]]) <- c("Rhat", "ESS")
    names(mcmc_diagnostics)[1:K] <- region_names
    names(mcmc_diagnostics)[K + 1] <- "Global"
  }else{
    # just get them for sigma2
    samples_matrix <- matrix(nrow = nstore, ncol = chains)
    for(chn in 1:chains){
      samples_matrix[,chn] <- MCMC_Output[[chn]][["V"]][3,] # first two rows are filler
    }
    mcmc_diagnostics[[K + 1]] <- matrix(convergence_diagnostics(samples_matrix),1,2)
    rownames(mcmc_diagnostics[[K + 1]]) <- vparam_names[3]
    colnames(mcmc_diagnostics[[K + 1]]) <- c("Rhat", "ESS")
    names(mcmc_diagnostics)[1:K] <- region_names
    names(mcmc_diagnostics)[K + 1] <- "Global"
  }
  print(mcmc_diagnostics)
  
  # Now extract and concatenate samples across chains
  if(K == 1){
    Xi <- matrix(nrow = chains*nstore, ncol = P)
    for(chn in 1:chains){
      start_idx_destination <- nstore*(chn-1) + 1
      end_idx_destination <- nstore*(chn)
      Xi[start_idx_destination:end_idx_destination,] <- t(MCMC_Output[[chn]][["Xi"]])
    }
    colnames(Xi) <- predictor_names
  }else{
    Xi <- array(dim = c(chains*nstore,P,K))
    for(k in 1:K){
      for(chn in 1:chains){
        start_idx_source <- P*(k - 1) + 1
        end_idx_source <- P*k
        start_idx_destination <- nstore*(chn-1) + 1
        end_idx_destination <- nstore*(chn)
        Xi[start_idx_destination:end_idx_destination,,k] <- t(MCMC_Output[[chn]][["Xi"]][start_idx_source:end_idx_source,])
      }
    }
    dimnames(Xi)[[2]] <- predictor_names
    dimnames(Xi)[[3]] <- region_names
  }

  # extract Xi0 (only if not single region style fit)
  if(!sr_style){
    for(chn in 1:chains){
      Xi0 <- matrix(nrow = chains*nstore, ncol = P)
      for(chn in 1:chains){
        start_idx_destination <- nstore*(chn-1) + 1
        end_idx_destination <- nstore*(chn)
        Xi0[start_idx_destination:end_idx_destination,] <- t(MCMC_Output[[chn]][["Xi0"]])
      }
    }
    colnames(Xi0) <- predictor_names
  }

  # extract Vparams
  if(!sr_style){
    for(chn in 1:chains){
      V <- matrix(nrow = chains*nstore, ncol = 3)
      for(chn in 1:chains){
        start_idx_destination <- nstore*(chn-1) + 1
        end_idx_destination <- nstore*(chn)
        V[start_idx_destination:end_idx_destination,] <- t(MCMC_Output[[chn]][["V"]])
      }
    }
    colnames(V) <- vparam_names
  }else{
    for(chn in 1:chains){
      V <- matrix(nrow = chains*nstore, ncol = 1)
      for(chn in 1:chains){
        start_idx_destination <- nstore*(chn-1) + 1
        end_idx_destination <- nstore*(chn)
        V[start_idx_destination:end_idx_destination] <- c(MCMC_Output[[chn]][["V"]][3,])
      }
    }
    colnames(V) <- "Variance(GP Random Effect)"
  }
  # Store samples in a list for output
  samples <- list()
  if(!sr_style){
    samples[["Xi"]] <- Xi
    samples[["Xi0"]] <- Xi0
    samples[["V"]] <- V
  }else{
    samples[["Xi"]] <- Xi
    samples[["V"]] <- V
  }
  
  summary_stats <- list()
  if(K == 1){
    summary_stats[[1]] = t(apply(samples$Xi, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x,c(0.025,0.975)))))
  }else{
    for(k in 1:K){
      summary_stats[[k]] = t(apply(samples$Xi[,,k], 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x,c(0.025,0.975)))))
    }
  }
  if(!sr_style){
    summary_stats[[K + 1]] = t(apply(samples$Xi0, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x,c(0.025,0.975)))))
    summary_stats[[K + 1]] = rbind(summary_stats[[K + 1]],t(apply(samples$V, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x,c(0.025,0.975))))))
  }else{
    for(k in 1:K){
      summary_stats[[k]] = rbind(summary_stats[[k]],t(apply(samples$V, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x,c(0.025,0.975))))))
    }
  }
  if(!sr_style){
    names(summary_stats) <- c(region_names,"Global")
  }else{
    names(summary_stats) <- region_names
  }
  print(summary_stats)
  
  # Construct a return object using 
  # - summary stats for each parameter
  # - rhat and ESS diagnostics for each parameter
  # - the list of samples compiled by smesir_mcmc
  # - summary stats for the "inferred number of cases" trajectory in each region
  # - summary stats for the "expected incidence event" trajectory in each region

  #beta_trajs <- design_matrices[[1]]%*%t(Xi)
  #matplot(t(apply(beta_trajs,1,function(samps) quantile(samps, c(0.025,0.975)))), type = "l", main = "Inferred Transmission Rate")
  output <- list(summary = summary_stats, 
                 mcmc_diagnostics = mcmc_diagnostics, 
                 samples = samples,
                 prior = prior,
                 formula = formula,
                 sr_style = sr_style,
                 design_matrices = design_matrices)
  return(output)
}