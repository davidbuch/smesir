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
#' @param formula Object of class "formula" 
#' @param data A named list containing the data with which the model will be fit. 
#' The list should include an entry for each term in the accompanying \code{formula} argument.
#' @param epi_params Epidemiologic parameters which must be specified by the user (typically obtained from side-information):
#' \code{region_populations} - Vector of populations of the regions studied, listed in the same order in which they are indexed in the data
#' \code{outbreak_times} - Vector of indices of the time interval at which the first cases are reported in each region
#' \code{mean_removal_time} - Average amount of time (in number or fractions of time intervals) that an individual remains infectious
#' \code{psi} - Vector of probabilities whose element \code{i} is the probability that the response event (case detection, death) occurs an infected individual \code{i - 1} time intervals after their infection
#' @param prior (Optional - reasonable default values are specified internally) A length 4 named list containing:
#' \code{ell} - lengthscale of the squared exponential kernal for the temporal random effect
#' \code{V0} - 3 variance hyperparameters for gaussian priors on the intercepts, coefficients, and temporal random effects;
#' \code{IGSR} - 3 pairs of shape and rate hyperparameters for inverse-gamma priors;
#' \code{expected_initial_infected_population} - the expected size of the infected population that appears at the beginning of the outbreak, used in an exponential prior; 
#' @param region_names Vector of names of the regions studied, listed in the same order in which they are indexed in the data
#' @export
smesir <- function(formula, data, epi_params, vaccinations = NULL, region_names = NULL, prior = NULL,
                   chains = 4, iter = 50000, warmup = 0, thin = max(floor((iter - warmup)/1000),1), 
                   min_adaptation_cycles = 5, min_samps_per_cycle = NULL, 
                   tempering_ratio = 0.2, quiet = TRUE, sr_style = NULL, 
                   seed = NULL){
  ## Unpack epi_params
  region_populations <- epi_params$region_populations
  outbreak_times <- epi_params$outbreak_times
  mean_removal_time <- epi_params$mean_removal_time
  if(length(dim(epi_params$incidence_probabilities)) != 2){
    K <- length(region_populations)
    epi_params$incidence_probabilities <- matrix(rep(epi_params$incidence_probabilities, K),ncol = K)
  }
  incidence_probabilities <- epi_params$incidence_probabilities
  
  #print(incidence_probabilities)
  
  if(!(is.null(vaccinations) || is.matrix(vaccinations))){
    stop("Vaccination data must be provided as a matrix 
         (rows = number of time points, cols = number of regions).")
  }
  #stop()
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
  if(K == 1){
    if(is.null(prior[["ell"]])) prior[["ell"]] <- J/5
    if(is.null(prior[["V0"]])) prior[["V0"]] <- c(10,10)
    if(is.null(prior[["expected_initial_infected_population"]])) prior[["expected_initial_infected_population"]] <- region_populations/10000
    if(is.null(prior[["IGSR"]])) prior[["IGSR"]] <- c(3,0.2)
    
    if(!is.numeric(prior[["ell"]]) || length(prior[["ell"]]) != 1 || prior[["ell"]] <= 0){
      stop("prior[['ell']] must be a positive scalar")
    }
    if(!is.numeric(prior[["V0"]]) || length(prior[["V0"]]) != 2 || any(prior[["ell"]] <= 0)){
      stop("prior[['V0']] must be a length 2 positive numeric vector (single-region style)")
    }
    if(!is.numeric(prior[["expected_initial_infected_population"]]) || length(prior[["expected_initial_infected_population"]]) != 1 || prior[["expected_initial_infected_population"]] <= 0){
      stop("prior[['expected_initial_infected_population']] must be a positive scalar")
    }
    if(!is.numeric(prior[["IGSR"]]) || length(prior[["IGSR"]]) != 2 || any(prior[["IGSR"]] <= 0)){
      stop("prior[['IGSR']] must be a length 2 positive numeric vector (single-region style)")
    }

    prior[["V0"]] <- c(prior[["V0"]],1) # last value placeholder
    prior[["IGSR"]] <- matrix(rep(prior[["IGSR"]],3), nrow = 3, ncol = 2,  byrow = TRUE)
  }else{
    if(is.null(prior[["ell"]])) prior[["ell"]] <- J/5
    if(is.null(prior[["V0"]])) prior[["V0"]] <- c(10,10,0.1)
    if(is.null(prior[["expected_initial_infected_population"]])) prior[["expected_initial_infected_population"]] <- mean(region_populations)/10000
    if(is.null(prior[["IGSR"]])) prior[["IGSR"]] <- matrix(c(rep(c(2.01,0.101),2),3,0.2), nrow = 3, ncol = 2, byrow = TRUE)
    
    if(!is.numeric(prior[["ell"]]) || length(prior[["ell"]]) != 1 || prior[["ell"]] <= 0){
      stop("prior[['ell']] must be a positive scalar")
    }
    if(!is.numeric(prior[["V0"]]) || length(prior[["V0"]]) != 3 || any(prior[["ell"]] <= 0)){
      stop("prior[['V0']] must be a length 3 positive numeric vector")
    }
    if(!is.numeric(prior[["expected_initial_infected_population"]]) || length(prior[["expected_initial_infected_population"]]) != 1 || prior[["expected_initial_infected_population"]] <= 0){
      stop("prior[['expected_initial_infected_population']] must be a positive scalar")
    }
    if(!is.numeric(prior[["IGSR"]]) || dim(prior[["IGSR"]]) != c(3,2) || any(prior[["IGSR"]] <= 0)){
      stop("prior[['IGSR']] must be a 3x2 positive numeric matrix")
    }
  }

  ell <- prior[["ell"]]
  V0 <- prior[["V0"]]
  IGSR <- prior[["IGSR"]]
  expected_initial_infected_population <- prior[["expected_initial_infected_population"]]
  
  if(is.null(vaccinations)){
    vaccinations <- matrix(0, nrow = J, ncol = K)
  }
  if(any(dim(vaccinations) != c(J,K)) || any(vaccinations > 1) || any(vaccinations < 0)){
    stop("'vaccinations' must be a JxK matrix of new vaccinations in each interval (as a percentage of each region's population).")
  }
  
  ## We will need the eigendecomposition of the scaled GP covariance matrix
  ## when we construct our design matrix
  
  # set 'r' here so it is the same for all regions
  S_ell <- exp(-as.matrix(dist(1:J, diag = TRUE, upper = TRUE)/ell)^2)
  kernel_decomp <- eigen(S_ell,symmetric = TRUE)
  r <- match(TRUE, cumsum(kernel_decomp$values)/sum(kernel_decomp$values) > .99)
  
  ## Construct Data Matrices for Model from user provided data list
  design_matrices <- list()
  tilde_off <- list() # matrices to adjust restricted regression parameter estimates
  response_matrix <- matrix(nrow = J, ncol = K)
  if(K == 1){
    covariate_dframe <- list()
    for(cname in covariate_names){
      if(cname %in% colnames(data)){
        covariate_dframe[[cname]] <- data[,cname]
      }else{
        stop(paste0("Covariate '",cname,"' not found in data."))
      }
    }
    covariate_dframe <- as.data.frame(covariate_dframe)
    design_mat <- model.matrix(as.formula(paste("~", as.character(formula)[3])), covariate_dframe)
    if(length(covariate_names) == 0){
      design_mat <- matrix(1, nrow = J, ncol = 1)
      colnames(design_mat) <- "(Intercept)"
    }
    nonblank <- which(apply(design_mat,2,function(a) any(a != 0)))
    X <- design_mat[,nonblank]
    pdm <- ncol(design_mat)
    TO <- diag(pdm+r) # adjustment matrix (tilde_off)
    TO[nonblank,(pdm+1):(pdm+r)] <- -solve(t(X)%*%X)%*%t(X)%*%kernel_decomp$vectors[,1:r]
    ImPx <- diag(nrow(X)) - X%*%solve(t(X)%*%X)%*%t(X) # Px_perp

    # restricted basis
    ebasis <- ImPx %*% kernel_decomp$vectors[,1:r]
    colnames(ebasis) <- paste("GP Basis Func.",1:r)
    vscales_theta <- kernel_decomp$values[1:r]
    
    tilde_off[[1]] <- TO
    design_matrices[[1]] <- cbind(design_mat,ebasis)
    response_matrix[,1] <- data[,response_name]
  }else{
    for(k in 1:K){
      covariate_dframe <- list()
      for(cname in covariate_names){
        if(cname %in% names(data)){
          covariate_dframe[[cname]] <- data[[cname]][,k]
        }else{
          stop(paste0("Covariate '",cname,"' column '",k,"' not found in data."))
        }
      }
      covariate_dframe <- as.data.frame(covariate_dframe)
      design_mat <- model.matrix(as.formula(paste( "~", as.character(formula)[3])), covariate_dframe)
      if(length(covariate_names) == 0){
        design_mat <- matrix(1, nrow = J, ncol = 1)
        colnames(design_mat) <- "(Intercept)"
      }
      nonblank <- which(apply(design_mat,2,function(a) any(a != 0)))
      X <- design_mat[,nonblank]
      pdm <- ncol(design_mat)
      TO <- diag(pdm+r) # adjustment matrix (tilde_off)
      TO[nonblank,(pdm+1):(pdm+r)] <- -solve(t(X)%*%X)%*%t(X)%*%kernel_decomp$vectors[,1:r] 
      ImPx <- diag(nrow(X)) - X%*%solve(t(X)%*%X)%*%t(X) # Px_perp

      # restricted basis
      ebasis <- ImPx %*% kernel_decomp$vectors[,1:r]
      colnames(ebasis) <- paste("GP Basis Func.",1:r)
      vscales_theta <- kernel_decomp$values[1:r]

      tilde_off[[k]] <-  TO
      design_matrices[[k]] <- cbind(design_mat,ebasis) ## attach eigenbasis
      response_matrix[,k] <- data[[response_name]][,k]
    }
  }
  if(anyNA(response_matrix) || anyNA(design_matrices, recursive = TRUE)){
    stop("Function 'smesir' is not currently set up to handle data with missing values.")
  }
  P <- ncol(design_matrices[[1]]) # number of "predictors" (including intercept and GP bases)
  nterms <- P - r - 1 # not the same as the number of covariates (e.g., factor expansions)
  predictor_names <- colnames(design_matrices[[1]])
  vparam_names <- c("Variance(Intercept)", "Variance(Covariate Coeffs.)", "Variance(GP Random Effect)")
  names(design_matrices) <- region_names
  colnames(response_matrix) <- region_names
  if(!sr_style){
    if(nterms == 0){
      V_used <- c(1,3)
    }else{
      V_used <- 1:3
    }
  }
  # when passing off to the MCMC functions, pass Y as a
  # matrix, whether or not K > 1. Similarly, pass XI as a matrix
  # and the Design_Matrices as a list (which is how I've already
  # constructed 'design_matrices' fortunately)
  if(!is.null(seed)) set.seed(seed)
  if(!quiet){
    print("Reached the Sampling Function")
  }
  if(is.null(min_samps_per_cycle)){
    min_samps_per_cycle <- 10*P*P # this should be pretty large since samples are autocorrelated
  }
  
  
  MCMC_Output <- smesir_mcmc(response_matrix, design_matrices, tilde_off, vaccinations, vscales_theta, V0, IGSR, 1/mean_removal_time, outbreak_times, expected_initial_infected_population,
                             region_populations, incidence_probabilities, min_adaptation_cycles, min_samps_per_cycle, chains,iter,warmup,thin,sr_style,quiet) # last arg is sr_style flag
  
  ## do convergence diagnostics here
  nstore <- floor((iter - warmup)/thin)
  mcmc_diagnostics <- list()
  for(k in 1:K){
    # get them for Xi
    mcmc_diagnostics[[k]] <- matrix(nrow = P + 2, ncol = 2)
    for(p in 1:P){
      samples_matrix <- matrix(nrow = nstore, ncol = chains)
      for(chn in 1:chains){
        source_row <- P*(k - 1) + p
        samples_matrix[,chn] <- MCMC_Output[[chn]][["Xi"]][source_row,]
      }
      mcmc_diagnostics[[k]][p,] <- convergence_diagnostics(samples_matrix)
    }
    # get them for IIP
    samples_matrix <- matrix(nrow = nstore, ncol = chains)
    for(chn in 1:chains){
      samples_matrix[,chn] <- MCMC_Output[[chn]][["IIP"]][k,]
    }
    mcmc_diagnostics[[k]][P + 1,] <- convergence_diagnostics(samples_matrix)
    # get them for DISP
    samples_matrix <- matrix(nrow = nstore, ncol = chains)
    for(chn in 1:chains){
      samples_matrix[,chn] <- MCMC_Output[[chn]][["DISP"]][k,]
    }
    mcmc_diagnostics[[k]][P + 2,] <- convergence_diagnostics(samples_matrix)
    
    rownames(mcmc_diagnostics[[k]]) <- c(predictor_names,"IIP","DISP")
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


  # Now extract and concatenate samples across chains
  if(K == 1){
    Xi <- matrix(nrow = chains*nstore, ncol = P)
    for(chn in 1:chains){
      start_idx_destination <- nstore*(chn-1) + 1
      end_idx_destination <- nstore*(chn)
      Xi[start_idx_destination:end_idx_destination,] <- t(MCMC_Output[[chn]][["Xi"]])
    }
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
    dimnames(Xi)[[3]] <- region_names
  }
  
  # Convert coefficient estimates and random effect basis back to original form
  # Transform back to unrestricted parameters and correct the basis
  # we will need the unaltered basis functions
  if(K == 1){
    Xi <- t(tilde_off[[1]] %*% t(Xi))
    design_matrices[[1]][,(P - r + 1):P] <- kernel_decomp$vectors[,1:r]
    colnames(Xi) <- predictor_names
  }else{
    for(k in 1:K){
      Xi[,,k] <- t(tilde_off[[k]] %*%t(Xi[,,k]))
      design_matrices[[k]][,(P - r + 1):P] <- kernel_decomp$vectors[,1:r]
    }
    dimnames(Xi)[[2]] <- predictor_names
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

  for(chn in 1:chains){
    IIP <- matrix(nrow = chains*nstore, ncol = K)
    for(chn in 1:chains){
      start_idx_destination <- nstore*(chn-1) + 1
      end_idx_destination <- nstore*(chn)
      IIP[start_idx_destination:end_idx_destination,] <- t(MCMC_Output[[chn]][["IIP"]])
    }
  }
  colnames(IIP) <- rep("IIP",K)
  
  for(chn in 1:chains){
    DISP <- matrix(nrow = chains*nstore, ncol = K)
    for(chn in 1:chains){
      start_idx_destination <- nstore*(chn-1) + 1
      end_idx_destination <- nstore*(chn)
      DISP[start_idx_destination:end_idx_destination,] <- t(MCMC_Output[[chn]][["DISP"]])
    }
  }
  colnames(DISP) <- rep("DISP",K)
  
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
    samples[["IIP"]] <- IIP
    samples[["DISP"]] <- DISP
    samples[["V"]] <- V
  }else{
    samples[["Xi"]] <- Xi
    samples[["IIP"]] <- IIP
    samples[["DISP"]] <- DISP
    samples[["V"]] <- V
  }
  
  # define a summarization function
  fournum <- function(x) c(mean = mean(x), sd = sd(x), quantile(x,c(0.025,0.975)))
  summary_stats <- list()
  if(K == 1){
    summary_stats[[1]] <- rbind(t(apply(matrix(samples$Xi[,1:(1+nterms)],ncol = 1+nterms), 2, fournum)),
                               fournum(samples$IIP),fournum(samples$DISP))
    rownames(summary_stats[[1]]) <- c(predictor_names[1:(1+nterms)],"IIP","DISP")
  }else{
    for(k in 1:K){
      summary_stats[[k]] <- rbind(t(apply(matrix(samples$Xi[,1:(1+nterms),k],ncol=1+nterms), 2, fournum)),
                                 fournum(samples$IIP[,k]),fournum(samples$DISP[,k]))
      rownames(summary_stats[[k]]) <- c(predictor_names[1:(1+nterms)],"IIP","DISP")
    }
  }
  if(!sr_style){
    summary_stats[[K + 1]] <- t(apply(matrix(samples$Xi0[,1:(1+nterms)],ncol=1+nterms), 2, fournum))
    rownames(summary_stats[[K + 1]]) <- predictor_names[1:(1+nterms)]
    summary_stats[[K + 1]] <- rbind(summary_stats[[K + 1]],t(apply(samples$V, 2, fournum)))
  }else{
    summary_stats[[K + 1]] <- t(apply(samples$V, 2, fournum))
  }
  if(!sr_style){
    names(summary_stats) <- c(region_names,"Global")
  }else{
    names(summary_stats) <- c(region_names,"Global")
  }
  if(!quiet){
    print(summary_stats)
  }
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
                 vaccinations = vaccinations,
                 prior = prior,
                 formula = formula,
                 sr_style = sr_style,
                 response_matrix = response_matrix,
                 design_matrices = design_matrices,
                 epi_params = epi_params,
                 region_names = region_names,
                 n_basis = r)
  return(output)
}