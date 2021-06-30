#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

inline double min(const double a, const double b){
  return (a + b - fabs(a - b))/2;
}
inline double max(const double a, const double b){
  return (a + b + fabs(a - b))/2;
}

arma::vec mvrnorm(arma::vec &mu, arma::mat &Sigma) {
  arma::vec X; X.randn(mu.size());
  arma::mat R = arma::chol(Sigma);
  return(trans(R) * X + mu);
}

//[[Rcpp::export]]
NumericVector solve_infections(const NumericVector beta, const double gamma, const int t_1, const double initial_impulse, const int N) {
  // Inputs:
  // t_1 - time of the outbreak
  // infection_impulse - compartment proportion discontinuity at the timestep before outbreak (o_k - 1)
  // beta - transmission rate values at timepoints o_k - 1, o_k, ..., t_J
  // gamma - inverse rate of removal (constant over time)
  // Output:
  // nu - sequence of new cases observed at timepoints o_k, o_k + 1, dots, t_J
  int J = beta.size();
  NumericVector nu(J); // will initialize with J zeros
  // first case detection occurs in the interval t_0 -> t_1
  int offset = t_1 - 1; // translate index to Cpp: nu_1 is labeled nu_0 here.
  nu(offset) = initial_impulse;
  double state[] = {1 - initial_impulse, initial_impulse};
  double infections = 0, removals = 0;
  for(int j = 1; j != (J - offset); ++j){
    infections = max(min(beta(offset + j)*state[0]*state[1],state[0]),0);
    removals = max(min(gamma*state[1], state[1]),0);
    state[0] = state[0] - infections;
    state[1] = state[1] + infections - removals;
    nu[offset + j] = infections;
  }
  nu[nu == 0] = 1e-16;
  return N*(nu + abs(nu))/2; // get positive part to remove round-off negativity 
}

//[[Rcpp::export]]
NumericVector solve_events(const NumericVector nu, // interval new infections 
                           const NumericVector psi // interval event probabilities
){
  int J = nu.size();
  int M = psi.size();
  int u = 0;
  NumericVector delta(J);
  for(int j = 0; j != J; ++j){
    u = ((j + 1) + M - abs(j + 1 - M))/2; // == min(j+1,M)
    for(int m = 0; m != u; ++m){
      delta[j] += nu[j - m] * psi[m];
    }
  }
  return delta;
}

//[[Rcpp::export]]
double log_poisd(const NumericVector y, const NumericVector lambda){
  double res = 0;
  for(int i = 0; i != y.length(); ++i){
    if(lambda[i] != 0 || y[i] != 0){ // if lambda[i] == 0, we have to be careful
      res += y[i]*log(lambda[i]) - lambda[i];
    }
  }
  return res;
}

//[[Rcpp::export]]
double log_llh(const NumericMatrix B, 
               const NumericMatrix Y, 
               const double gamma, // inverse infectious period
               const IntegerVector T_1, // outbreak indices 
               const NumericVector Initial_Impulse, // initial impulses of infection
               const IntegerVector N, // region populations
               const NumericVector psi, // interval event probability
               const NumericMatrix frailty
){
  const int K = Y.ncol();//, J = Y.nrow();
  double llh = 0;
  for(int k = 0; k != K; ++k){
    llh += log_poisd(Y(_,k),
                     frailty(_,k)*solve_events(
                       solve_infections(B(_,k),gamma,T_1(k),Initial_Impulse(k),N(k)),
                       psi)
    );
  }
  return llh;
}

double log_prior(const arma::mat & Xi,   // PxK matrix of local parameters
                 const arma::vec & Xi0,     // P vector of global parameters
                 const arma::vec & V,      // P? vector of variances
                 const arma::mat & B
){
  // check for B out of bounds before computing the rest of the prior
  bool B_out_of_bounds = arma::any(arma::vectorise(B) < 0);
  if(B_out_of_bounds){
    return(log(0));
  }
  arma::mat X(Xi.each_col() - Xi0);
  // assume inverse variance diagonal is fixed in advance
  // we actually don't need to add the prior contributions
  // of the global parameters or the variance terms because 
  // we use update those parameters with gibbs proposals
  return( -0.5 * arma::as_scalar( arma::dot(((X*(X.t())).eval()).diag(), 1/V) ) );
}

double log_posterior(const arma::mat & Xi,      // PxK matrix of local parameters
                     const arma::vec & Xi0,     // P vector of global parameters
                     const arma::vec & V,       // P? vector of variances
                     const NumericMatrix Y, // auto pass by reference, I believe
                     const List Design_Matrices,  // List of K design matrices
                     const double gamma, // inverse infectious period
                     const IntegerVector T_1, // outbreak indices 
                     const NumericVector Initial_Impulse, // initial impulses of infection
                     const IntegerVector N, // region populations
                     const NumericVector psi, // interval event probability
                     const NumericMatrix frailty
)
{
  const int J = Y.nrow(), K = Y.ncol(); //, P = Xi.n_rows;
  arma::mat B(J,K);
  for(int k = 0; k != K; ++k){
    arma::mat dmat = Design_Matrices[k];
    B.col(k) = dmat * Xi.col(k);
  }
  double res = 0;
  res += log_llh(NumericMatrix(J,K,B.begin()), Y, gamma, T_1, Initial_Impulse, N, psi, frailty);
  //Rcout << res << " log llh \n" ;
  res += log_prior(Xi,Xi0,V,B);
  //Rcout << res << " log llh + log prior \n";
  return(res);
}

void update_frailty(NumericMatrix frailty,
                      const double dispersion,
                      const NumericMatrix Y,
                      const arma::mat & Xi,
                      const List Design_Matrices,
                      const double gamma,
                      const IntegerVector T_1,
                      const NumericVector Initial_Impulse,
                      const IntegerVector N,
                      const NumericVector psi)
{
  const int J = Y.nrow(), K = Y.ncol(); //, P = Xi.n_rows;
  arma::mat B(J,K);
  for(int k = 0; k != K; ++k){
    arma::mat dmat = Design_Matrices[k];
    B.col(k) = dmat * Xi.col(k);
  }
  NumericMatrix expected_events(J,K);
  for(int k = 0; k != K; ++k){
    expected_events(_,k) = solve_events(solve_infections(as<NumericVector>(wrap(B.col(k))),gamma,T_1(k),Initial_Impulse(k),N(k)),psi);
  }
  
  for(int j = 0; j != J; ++j){
    for(int k = 0; k != K; ++k){
      frailty(j,k) = R::rgamma((Y(j,k) + 1/dispersion),1/(expected_events(j,k) + 1/dispersion));
    }
  }
}

void update_xi0(arma::mat const & Xi,
                arma::vec & Xi0,
                arma::vec const & V,
                arma::vec const & V0
)
{
  int K = Xi.n_cols;
  arma::mat post_var = arma::diagmat(1/(K*(1/V) + (1/V0))); // element-wise inverse
  arma::vec post_mean = post_var * (K/V) % mean(Xi,1); // rowMeans(Xi)
  Xi0 = mvrnorm(post_mean, post_var);
}

void update_Vparam(arma::mat const & Xi,
                   arma::vec const & Xi0,
                   arma::vec & Vparam,
                   NumericMatrix IGSR, // Inverse-Gamma Shape and Rate hyperparams
                   arma::vec const & lambda,
                   int const ncovs,
                   int const nbases,
                   bool const vparam_key[]
)
{
  arma::mat X(Xi.each_col() - Xi0);
  int K = X.n_cols;
  int P = X.n_rows; // 1 + ncovs + nbases
  arma::vec ss_diag(P);
  for(int p = 0; p != P; ++p){
    ss_diag(p) = arma::dot(X.row(p),X.row(p));
  }
  
  // update regional intercept variance
  double a = 0, b = 0;
  if(vparam_key[0]){
    a = (2*IGSR(0,0) + K)/2;
    b = (2*IGSR(0,1) + ss_diag(0))/2;
    Vparam(0) = 1/arma::as_scalar(arma::randg<arma::vec>(1,arma::distr_param(a,b)));
  }
  // update regional covariate coeff variance
  if(vparam_key[1]){
    a = (2*IGSR(1,0) + ncovs*K)/2;
    /* subvec indices are inclusive */
    b = (2*IGSR(1,1) + arma::sum(ss_diag.subvec(1,1+ncovs-1)) )/2;
    Vparam(1) = 1/arma::as_scalar(arma::randg<arma::vec>(1,arma::distr_param(a,b)));
  }
  
  //update regional gp basis coeff variance
  if(vparam_key[2]){
    a = (2*IGSR(2,0) + nbases*K)/2;
    /* subvec indices are inclusive */
    b = (2*IGSR(2,1) + arma::sum(ss_diag.subvec(1+ncovs,P-1)/lambda))/2;
    Vparam(2) = 1/arma::as_scalar(arma::randg<arma::vec>(1,arma::distr_param(a,b)));
  }
}

void expand_Vparam(arma::vec & V,
                   arma::vec const & Vparam,
                   arma::vec const & lambda,
                   int const ncovs,
                   int const nbases
)
{
  V(0) = Vparam(0); // intercept prior variance
  for(int i = 1; i < 1 + ncovs; ++i){
    V(i) = Vparam(1); // covariate coefficient prior variance
  }
  for(int i = 1 + ncovs; i < 1 + ncovs + nbases; ++i){
    V(i) = Vparam(2)*lambda(i - (1 + ncovs)); // covariate coefficient prior variance
  }
}


//[[Rcpp::export]]
List smesir_mcmc(const NumericMatrix Y,
                 const List Design_Matrices,  // List of K design matrices
                 const arma::vec & lambda, // GP basis eigenvalues
                 const arma::vec & V0param, // pre-expanded variance hyperparameters
                 const NumericMatrix IGSR, // (3x2 if hierarchical, 1x2 ow) Gamma Shape and Rate hyperparameters
                 const double gamma, // inverse infectious period
                 const IntegerVector T_1, // outbreak indices 
                 const NumericVector Initial_Impulse, // initial impulses of infection
                 const double dispersion, // dispersion == 0 implies poisson likelihood
                 const IntegerVector N, // region populations
                 const NumericVector psi, // interval event probability
                 const double tempering_ratio,
                 const int ncycles,
                 const int samps_per_cycle,
                 const int nchain, 
                 const int iter, 
                 const int warmup, 
                 const int thin,
                 const bool sr_style,
                 const bool quiet
)
{
  if(!quiet){
    Rcout << "Reached the inside of the function... \n";
  }
  const int K = Y.ncol(), J = Y.nrow();
  NumericMatrix const & Eg_Design_Mat = Design_Matrices[0];
  const int P = Eg_Design_Mat.ncol();
  
  const int nbases = lambda.n_elem; //number of basis coefficients
  const int ncovs = P - (nbases + 1); // number of covariate coefficients
  const bool vparam_key[] = {!sr_style, (ncovs > 0) && !sr_style, nbases > 0};
  
  arma::vec V0(P); 
  expand_Vparam(V0,V0param,lambda,ncovs,nbases);
  
  const int nstore = (iter - warmup)/thin; // automatically floored
  
  if(!quiet){
    Rcout << "Begin iterating over chains... \n";
  }
  
  List MCMC_output(nchain);
  for(int chn = 0; chn != nchain; ++chn){
    if(!quiet){
      Rcout << "Chain " << chn + 1 << " of " << nchain <<":\n";
      Rcout << "Initializing local coefficients... \n\tRegions:";
    }
    // randomly initialize parameters
    arma::mat Xi(P,K);
    arma::vec Xik;
    for(int k = 0; k != K; ++k){
      arma::mat dmat = Design_Matrices[k];
      bool within_prior_support = false;
      int counter = 0;
      while(!within_prior_support){
        Xik.randn(P);
        within_prior_support = arma::all(dmat*Xik > 0);
        counter = (1 + counter) % 1000;
        if(counter == 0){checkUserInterrupt();}
      }
      Xi.col(k) = Xik;
      if(!quiet){
        Rcout << " " << k + 1;
      }
    }

    NumericMatrix frailty(J,K); // initialize from prior
    for(int j = 0; j != J; ++j){
      for(int k = 0; k != K; ++k){
        frailty(j,k) = R::rgamma(1/dispersion,dispersion);
      }
    }
    
    if(!quiet){
      Rcout << "\n...done!\n";
      Rcout << "Initializing global coefficients and variance parameters...";
    }
    arma::vec V(P);
    arma::vec Vparam(3); // initialize from prior distribution
    for(int i = 0; i != 3; ++i){
      if(vparam_key[i]){
        Vparam(i) = 1/arma::as_scalar(arma::randg<arma::vec>(1,arma::distr_param(IGSR(i,0),IGSR(i,1))));
      }else{
        Vparam(i) = V0param(i);
      }
    }
    expand_Vparam(V,Vparam,lambda,ncovs,nbases);
    
    // set Xi0 to vec(0) and keep it that way if sr_style = true
    arma::vec Xi0(P,arma::fill::zeros); 
    if(!sr_style){
      update_xi0(Xi,Xi0,V,V0); // go ahead and initialize Xi0 among the Xi
    }
    arma::mat Xi_prop = Xi; // armadillo does a deep copy
 
    if(!quiet){
      Rcout << " done!\n";
      // Rcout << "\t\t Initialized Xi0 = \n";
      // Rcout << "\t\t " << Xi0.t() << "\n";
      // Rcout << "\t\t Initialized (tau_alpha, tau_phi, sigma^2) = \n";
      // Rcout << "\t\t" << Vparam.t() << "\n";
      Rcout << "Begin adaptation... \n\tAdaptation cycle:";
    }
    
    /* Adaptation Phase */
    arma::vec C0diag(P,arma::fill::ones);
    C0diag.subvec(1+ncovs,P - 1) = lambda/10; C0diag = (0.1/P)*C0diag;
    arma::mat C0 = arma::diagmat(C0diag);
    arma::cube R(P,P,K);
    for(int k = 0; k != K; ++k){
      R.slice(k) = arma::chol(C0);
    }
    // /*
    // arguments tempering_ratio, ncycles, samps per cycle, quiet
    IntegerVector xi_samp_counts(K,0);
    IntegerVector xi_samps_since_last_accepted(K,0);
    NumericVector acc_rates(K,0.0);
    int adaptation_cycle = 1;
    bool cycle_complete = false;
    double tempering_factor = std::pow(tempering_ratio,ncycles - 1);
    double tempering_increment = 1.0; // this will be set later based on the tempering factor at the end of the first cycle
    bool adapting = true;
    bool achieved_optimal_acc_rate = false;
    int it = 0; // within-adaptation iteration count
    int cycles_post_tempering = 0; // count cycles after tempering_factor == 1
    arma::cube ap_Xi(samps_per_cycle,P,K);
    while(adapting){
      it += 1; // update iteration number
      if((it % 1000) == 0){checkUserInterrupt();}
      for(int k = 0; k != K; ++k){
        if(xi_samp_counts[k] < samps_per_cycle){
          Xik = Xi.col(k);
          Xi_prop.col(k) = Xik + arma::trans(R.slice(k)) * arma::randn(P);
          double log_acc_prob = log_posterior(Xi_prop,Xi0,V,Y,Design_Matrices,gamma,T_1,Initial_Impulse,N,psi,frailty) - 
            log_posterior(Xi,Xi0,V,Y,Design_Matrices,gamma,T_1,Initial_Impulse,N,psi,frailty);
          if(arma::as_scalar(arma::randu(1)) < std::exp(log_acc_prob)){
            Xi.col(k) = Xi_prop.col(k);
            ap_Xi(xi_samp_counts[k],0,k,arma::size(1,P,1)) = Xi_prop.col(k);
            xi_samp_counts[k] += 1;
            xi_samps_since_last_accepted[k] = 0;
            acc_rates[k] = (1 + (it - 1)*acc_rates[k])/it;
          }else{
            if(tempering_factor == 1){ // stop being greedy
              ap_Xi(xi_samp_counts[k],0,k,arma::size(1,P,1)) = Xi.col(k);
              xi_samp_counts[k] += 1;
            }
            xi_samps_since_last_accepted[k] += 1;
            acc_rates[k] = (0 + (it - 1)*acc_rates[k])/it;
            Xi_prop.col(k) = Xi.col(k);
          }
          if(xi_samps_since_last_accepted[k] > 100){
            if(adaptation_cycle == 1){
              tempering_factor = tempering_factor/5;
            }
            R.slice(k) = R.slice(k)/3.16; // divide proposal covariance by 10
            for(int kk = 0; kk != K; ++kk){
              xi_samps_since_last_accepted[kk] = 0;
            }
          }
        }
      }
      
      // global params with gibbs proposals
      if(!sr_style){
        update_xi0(Xi,Xi0,V,V0); // modifies Xi0 in place
      }

      update_Vparam(Xi,Xi0,Vparam,IGSR,lambda,ncovs,nbases,vparam_key);
      expand_Vparam(V,Vparam,lambda,ncovs,nbases);

      update_frailty(frailty,dispersion,Y,Xi,Design_Matrices,gamma,T_1,Initial_Impulse,N,psi);

      cycle_complete = true;
      for(int k = 0; k!=K; ++k){
        cycle_complete = cycle_complete && (xi_samp_counts[k] == samps_per_cycle);
      }
      if(cycle_complete){
        //compute covariances
        for(int k = 0; k!= K; ++k){
          // magic number comes from 2.38^2, see article Haario et al (2001)
          if(tempering_factor < 1){
            R.slice(k) = arma::chol(5.66*arma::cov(ap_Xi.slice(k))/P);
          }else{ // after tempering, begin accumulating samples
            arma::mat oldSigHat = arma::trans(R.slice(k)) * R.slice(k);
            arma::mat newSigHat = 5.66*arma::cov(ap_Xi.slice(k))/P;
            arma::mat meanSigHat = ((cycles_post_tempering*samps_per_cycle - P)*oldSigHat + 
              (samps_per_cycle - P)*newSigHat)/((cycles_post_tempering+1)*samps_per_cycle - P);
            R.slice(k) = arma::chol(meanSigHat);
            cycles_post_tempering += 1;
          }
          
        }
        
        // print status updates
        if(!quiet){
          Rcout << " " << adaptation_cycle;
          Rcout << "Adaptation Cycle No.: " << adaptation_cycle << "\n";
          Rcout << "Total Iterations: " << it << "\n";
          Rcout << "Acceptance Rates: " << acc_rates << "\n";
          Rcout << "Tempering Factor: " << tempering_factor << "\n";
        }          
        
        // was this the last cycle?
        achieved_optimal_acc_rate = true;
        for(int k = 0; k!= K; ++k){
          achieved_optimal_acc_rate = achieved_optimal_acc_rate && (0.25 < acc_rates[k]) && (acc_rates[k] < 0.4);
          if(acc_rates[k] > 0.6 && tempering_factor == 1){
            R.slice(k) = 3.16*R.slice(k); // multiply covariance by 10
          }
        }
        
        if((adaptation_cycle > ncycles) && ((cycles_post_tempering >= 10) || achieved_optimal_acc_rate)){
          adapting = 0;
        }else{
          it = 0; // reset iteration number
          if(adaptation_cycle == 1){
            tempering_increment = (1 - tempering_factor)/ncycles;
          }
          adaptation_cycle += 1;
          tempering_factor += tempering_increment;
          tempering_factor = std::min(tempering_factor, 1.0);
          for(int k = 0; k != K; ++k){
            xi_samp_counts[k] = 0;
            acc_rates[k] = 0.0;
          }
        }
        if(adaptation_cycle > 100){
          stop("Over 100 adaptation cycles required: is something wrong?");
        }
      }
    }
    if(!quiet){
      Rcout << "\n...done!\n";
    }
    // */

    /* Main Sampling Phase */
    
    // allocate storage
    arma::mat chain_Xi(P*K,nstore);
    arma::mat chain_Xi0(P,nstore);
    arma::mat chain_V(3,nstore);
    
    Rcout << "Begin Sampling... \n\tPercent complete:";
    int percent_complete = 0;
    int percent_increment = 10;
    double sample_increment = iter/percent_increment;
    double completion_checkpoint = sample_increment;
    for(int it = 0; it != iter; ++it){
      if((it % 1000) == 0){checkUserInterrupt();}
      if(it > completion_checkpoint){
        completion_checkpoint += sample_increment;
        percent_complete += percent_increment;
        Rcout << " " << percent_complete;
      }
      // update local params
      // To-Do: Write a matrix-style update and profile to compare
      for(int k = 0; k != K; ++k){
        Xik = Xi.col(k);
        Xi_prop.col(k) = Xik + arma::trans(R.slice(k)) * arma::randn(P);
        double log_acc_prob = log_posterior(Xi_prop,Xi0,V,Y,Design_Matrices,gamma,T_1,Initial_Impulse,N,psi,frailty) - 
          log_posterior(Xi,Xi0,V,Y,Design_Matrices,gamma,T_1,Initial_Impulse,N,psi,frailty);
        if(arma::as_scalar(arma::randu(1)) < std::exp(log_acc_prob)){
          Xi.col(k) = Xi_prop.col(k);
        }else{
          Xi_prop.col(k) = Xi.col(k);
        }
      }
      
      update_frailty(frailty,dispersion,Y,Xi,Design_Matrices,gamma,T_1,Initial_Impulse,N,psi);
      
      // global params with gibbs proposals
      if(!sr_style){
        update_xi0(Xi,Xi0,V,V0); // modifies Xi0 in place
      }
      
      update_Vparam(Xi,Xi0,Vparam,IGSR,lambda,ncovs,nbases,vparam_key);
      expand_Vparam(V,Vparam,lambda,ncovs,nbases);
      
      if(!(it < warmup) && !(it % thin)){
        int s = (it - warmup)/thin;
        chain_Xi.col(s) = Xi.as_col();
        chain_Xi0.col(s) = Xi0;
        chain_V.col(s) = Vparam;
      }
    }
    Rcout << " 100\n...done!\n\n";
    List chain_output;
    chain_output["Xi"] = chain_Xi;
    chain_output["Xi0"] = chain_Xi0;
    chain_output["V"] = chain_V;
    MCMC_output[chn] = chain_output;
  }
  return(MCMC_output);
}

// If we use Rcpp::sourceCpp
// this R code will be automatically run after the compilation.

/*** R
# # Simulate Data
# J <- 50; K <- 5
# T_1 <- c(5,1,3,3,1) # first case reported at the end of interval t_1.
# N <- c(1e5,1e6,5e5,1e4,5e5)
# Initial_Impulse <- 10/N
# P <- 2
# twomat <- function(d,e){ h <- c(rep(0,d),rep(1,e - d)); return(matrix(c(rep(1,e),h),ncol = 2))}
# Design_Matrices <- list(twomat(8,J),twomat(3,J),twomat(5,J),twomat(10,J),twomat(3,J))
# Xi0 <- c(3,-1.5)
# Xi <- matrix(rnorm(P*K,rep(Xi0,K),sd = 0.2), ncol  = K)
# gamma <- 2/3
# psi <- 0.01*c(0.25,0.5,0.25)
# lambda <- 1
# V0 <- c(10,10,1) # sigma^2_0 should be small
# IGSR <- matrix(rep(2,6),nrow = 3, ncol = 2)
# Y <- matrix(nrow = J, ncol = K)
# for(k in 1:K){
#   Y[,k] <- rpois(J,solve_events(solve_infections(Design_Matrices[[k]]%*%Xi[,k],
#                                                  gamma, T_1[k],
#                                                  Initial_Impulse[k], N[k]),psi))
# 
# }
# 
# matplot(Y, type = "l", xlab = "day", ylab = NA, main = "Simulated Deaths")
# 
# MCMC_Output <- smesir_mcmc(Y,Design_Matrices,lambda,V0,IGSR,gamma,T_1, Initial_Impulse,
#                            N,psi,0.2,5,10*P*P,2,1000,500,5,FALSE,FALSE) # last arg is sr_style flag
# 
# matplot(t(MCMC_Output[[1]][["Xi"]][c(1,3,5,7,9),]),type = "l", ylab = NA)
# matplot(t(MCMC_Output[[1]][["Xi"]][c(2,4,6,8,10),]),type = "l", ylab = NA)
# trajs <- Y
# for(k in 1:K){
#   trajs <- cbind(trajs,solve_events(solve_infections(
#     Design_Matrices[[k]]%*%rowMeans(MCMC_Output[[1]][["Xi"]][(2*(k - 1) + 1):(2*(k  - 1) + 2),]), gamma, T_1[k], Initial_Impulse[k], N[k]),psi))
# }
# matplot(trajs, type = "l",ylab = NA)
# summary(apply(t(MCMC_Output[[1]][["V"]]),2,sqrt))


*/
