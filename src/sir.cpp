#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

inline double min(const double a, const double b){
  return (a + b - fabs(a - b))/2;
}
inline double max(const double a, const double b){
  return (a + b + fabs(a - b))/2;
}

double log_poisd(const NumericVector y, const NumericVector lambda){
  double res = 0;
  for(int i = 0; i != y.length(); ++i){
    if(lambda[i] != 0 || y[i] != 0){ // if lambda[i] == 0, we have to be careful
      res += y[i]*log(lambda[i]) - lambda[i];
    }
  }
  return res;
}

arma::vec mvrnorm(arma::vec &mu, arma::mat &Sigma) {
  arma::vec X; X.randn(mu.size());
  arma::mat R = arma::chol(Sigma);
  return(trans(R) * X + mu);
}

//[[Rcpp::export]]
NumericVector solve_infections(const NumericVector beta, // time-varying transmission rate
                               const double gamma, // inverse of average infectious period
                               const int T_1,   // time index of region's outbreak
                               const double II, // region's initial impulse of infection
                               const int N // region population
){
  int J = beta.size();
  NumericVector nu(J); // initialize to length J vector of 0s
  
  // first case detection occurs in the interval t_0 -> t_1
  int offset = T_1 - 1; // translate index to Cpp: nu_1 is labeled nu_0 here.
  nu(offset) = II;
  double state[] = {1 - II, II};
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


// single region's contribution to the log posterior density
double log_posterior_single(const arma::vec & Xi, // P vector of region's parameters
                     const arma::vec & Xi0,     // P vector of global parameters
                     const arma::vec & V,       // variance parameters
                     const NumericVector Y, // region's event time series
                     const arma::mat dmat, // region's design matrices
                     const double gamma, // inverse average infectious period
                     const int T_1, // region's outbreak index
                     const double II, // region's initial impulse
                     const int N, // region population
                     const NumericVector psi // interval event probability
){
  arma::vec B = dmat * Xi;
  // check if parameters are within prior support first
  bool B_out_of_bounds = arma::any(arma::vectorise(B) < 0);
  if(B_out_of_bounds){
    return(log(0));
  }
  
  double res = 0;
  // add log prior contribution
  arma::mat X(Xi - Xi0);
  res += -0.5 * arma::as_scalar( arma::dot(((X*(X.t())).eval()).diag(), 1/V) );
  
  // add log likelihood contribution
  res += log_poisd(Y,solve_events(solve_infections(as<NumericVector>(wrap(B)),gamma,T_1,II,N),psi));

  return(res);
}

void update_ii(const arma::mat & Xi,      // PxK matrix of local parameters
               const NumericMatrix Y, // auto pass by reference, I believe
               const List Design_Matrices,  // List of K design matrices
               const double gamma, // inverse infectious period
               const IntegerVector T_1, // outbreak indices 
               arma::vec & II, // initial impulses of infection
               const IntegerVector N, // region populations
               const NumericVector psi, // interval event probability
               const double prior_rate_param,
               const arma::vec & proposal_sd,
               const double tempering_factor
               )
{
  const int K = Y.ncol();
  
  for(int k = 0; k != K; ++k){
    arma::mat dmat = Design_Matrices[k];
    arma::vec B = dmat * Xi.col(k);
    
    double prop = II[k] + R::rnorm(0,proposal_sd[k])/N[k];
    prop = (prop > 0 ? prop : -prop); // reflect at 0
    double ac_prob = std::exp(
      - (prop - II[k]) * (N[k] * prior_rate_param) + // log prior ratio plus log likelihood ratio
      tempering_factor * (log_poisd(Y(_,k),solve_events(solve_infections(as<NumericVector>(wrap(B)),gamma,T_1[k],prop,N[k]),psi)) - 
        log_poisd(Y(_,k),solve_events(solve_infections(as<NumericVector>(wrap(B)),gamma,T_1[k],II[k],N[k]),psi))));
    
    if(R::runif(0,1) < ac_prob){
      II[k] = prop;
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
  
  // update variance of local intercept params around the global intercept param
  double a = 0, b = 0;
  if(vparam_key[0]){
    a = (2*IGSR(0,0) + K)/2;
    b = 2/(2*IGSR(0,1) + ss_diag(0)); // scale parameter
    Vparam(0) = 1/arma::randg<double>(arma::distr_param(a,b));
  }
  // update variance of local coefficient params around the global coefficient params
  if(vparam_key[1]){
    a = (2*IGSR(1,0) + ncovs*K)/2;
    /* subvec indices are inclusive */
    b = 2/(2*IGSR(1,1) + arma::sum(ss_diag.subvec(1,1+ncovs-1))); // scale parameter
    Vparam(1) = 1/arma::randg<double>(arma::distr_param(a,b));
  }
  
  //update variance of local gp basis params around global gp basis params
  if(vparam_key[2]){
    a = (2*IGSR(2,0) + nbases*K)/2;
    /* subvec indices are inclusive */
    b = 2/(2*IGSR(2,1) + arma::sum(ss_diag.subvec(1+ncovs,P-1)/lambda)); // scale parameter
    Vparam(2) = 1/arma::randg<double>(arma::distr_param(a,b));
  }
}

// expand variance params to allow vectorized operations in the log posterior function
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
                 const double prior_expected_ip, // expected infectious population
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
  const int K = Y.ncol(); //, J = Y.nrow();
  
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
    arma::vec Xik_prop = Xik; // armadillo does a deep copy
    

    // initialize initial impulse
    const double ii_prior_rate = 1.0/prior_expected_ip;
    arma::vec ii_proposal_sd(K,arma::fill::ones); // we will adapt this based on sample variance later
    arma::vec II(K);
    for(int k = 0; k != K; ++k){
      II[k] = R::rgamma(1,1/(ii_prior_rate*N[k])); // randomly initialize II ~ exp(1/(pr*N[k]))
    }
    
    if(!quiet){
      Rcout << "\n...done!\n";
      Rcout << "Initializing global coefficients and variance parameters...";
    }
    arma::vec V(P);
    arma::vec Vparam(3); // initialize from prior distribution
    for(int i = 0; i != 3; ++i){
      if(vparam_key[i]){
        Vparam(i) = 1/arma::randg<double>(arma::distr_param(IGSR(i,0),1/IGSR(i,1))); // shape, scale parameterization
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

    if(!quiet){
      Rcout << " done!\n";
      Rcout << "Begin adaptation... \n\tAdaptation cycle:";
    }
    
    /* Adaptation Phase */
    
    // prep to adapt Xi proposal covariance
    arma::vec C0diag(P,arma::fill::ones);
    C0diag.subvec(1+ncovs,P - 1) = lambda/10; C0diag = (0.1/P)*C0diag;
    arma::mat C0 = arma::diagmat(C0diag);
    arma::cube R(P,P,K);
    for(int k = 0; k != K; ++k){
      R.slice(k) = arma::chol(C0);
    }
    
    // prep to adapt II proposal variances
    arma::vec ii_current_mean(K,arma::fill::ones);
    arma::vec ii_current_var(K,arma::fill::ones);
    int ii_samp_count = 0;
    
    // arguments tempering_ratio, ncycles, samps per cycle, quiet
    IntegerVector xi_samp_counts(K,0);
    IntegerVector xi_samps_since_last_accepted(K,0);
    NumericVector acc_rates(K,0.0);
    int completed_cycle_count = 0;
    double tempering_factor = (double) completed_cycle_count/ncycles;
    bool adapting = true;
    int it = 0; // within-cycle iteration count
    arma::cube ap_Xi(samps_per_cycle,P,K);
    while(adapting){
      it += 1; // update iteration number
      if((it % 1000) == 0){checkUserInterrupt();}
      for(int k = 0; k != K; ++k){
        if(xi_samp_counts[k] < samps_per_cycle){
          Xik = Xi.col(k);
          Xik_prop = Xik + arma::trans(R.slice(k)) * arma::randn(P);
          double log_acc_prob = log_posterior_single(Xik_prop,Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),II(k),N(k),psi) - 
            log_posterior_single(Xik,Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),II(k),N(k),psi);
          if(R::runif(0,1) < std::exp(log_acc_prob)){
            Xi.col(k) = Xik_prop;
            ap_Xi(xi_samp_counts[k],0,k,arma::size(1,P,1)) = Xik_prop;
            xi_samp_counts[k] += 1;
            xi_samps_since_last_accepted[k] = 0;
            acc_rates[k] = (1 + (it - 1)*acc_rates[k])/it;
          }else{
            if(completed_cycle_count > ncycles){ // stop being greedy
              ap_Xi(xi_samp_counts[k],0,k,arma::size(1,P,1)) = Xi.col(k);
              xi_samp_counts[k] += 1;
            }
            xi_samps_since_last_accepted[k] += 1;
            acc_rates[k] = (0 + (it - 1)*acc_rates[k])/it;
          }
          if(xi_samps_since_last_accepted[k] > 100){
            R.slice(k) = R.slice(k)/3.16; // divide proposal covariance by 10
            xi_samps_since_last_accepted[k] = 0;
          }
        }
      }
      
      // global params with gibbs proposals
      if(!sr_style){
        update_xi0(Xi,Xi0,V,V0); // modifies Xi0 in place
      }
      
      // update II
      update_ii(Xi,Y,Design_Matrices,gamma,T_1,II,N,psi,
                ii_prior_rate,ii_proposal_sd,tempering_factor);
      // update II proposal standard deviation
      ii_samp_count += 1;
      for(int k = 0; k != K; ++k){ 
        ii_current_var[k] = ((ii_samp_count - 1)*ii_current_var[k]/ii_samp_count) + 
          ii_current_mean[k]*ii_current_mean[k] + N[k]*II[k]*N[k]*II[k]/ii_samp_count; // partially update var
        ii_current_mean[k] = (ii_samp_count*ii_current_mean[k] + N[k]*II[k])/(ii_samp_count + 1); // update mean
        ii_current_var[k] -= (ii_samp_count + 1)*ii_current_mean[k]*ii_current_mean[k]/ii_samp_count; // finish updating var
        ii_proposal_sd[k] = 2.4*std::sqrt(ii_current_var[k]);
      }
      
      // update variance parameters
      update_Vparam(Xi,Xi0,Vparam,IGSR,lambda,ncovs,nbases,vparam_key);
      expand_Vparam(V,Vparam,lambda,ncovs,nbases);

      // is the adaptation cycle complete?
      if(is_true(all(xi_samp_counts == samps_per_cycle))){
        // increment the completed cycle count
        completed_cycle_count += 1;
        
        // do we need to keep adapting? (are we short on cycles, or are acc_rates unacceptable)
        adapting = (completed_cycle_count < ncycles) || is_true(any(acc_rates < 0.1)) || is_true(any(acc_rates > 0.4));
        
        // compute new Xi proposal covariances
        for(int k = 0; k!= K; ++k){
          // magic number comes from 2.38^2, see article Haario et al (2001)
          if(completed_cycle_count < ncycles){
            R.slice(k) = arma::chol(5.66*arma::cov(ap_Xi.slice(k))/P);
          }else{ // after ncycle, begin retaining covariance information across cycles
            int extra_cycle_count = completed_cycle_count - (ncycles - 1);
            arma::mat oldSigHat = arma::trans(R.slice(k)) * R.slice(k);
            arma::mat newSigHat = 5.66*arma::cov(ap_Xi.slice(k))/P;
            arma::mat meanSigHat = ((extra_cycle_count*samps_per_cycle - P)*oldSigHat + 
              (samps_per_cycle - P)*newSigHat)/((extra_cycle_count+1)*samps_per_cycle - P);
            R.slice(k) = arma::chol(meanSigHat);
          }
        }
        // update tempering factor for II
        tempering_factor = std::min((double) completed_cycle_count/ncycles, 1.0);
        
        // print status updates
        if(!quiet){
          Rcout << "Adaptation Cycle: " << completed_cycle_count << "\n";
          Rcout << "Total Iterations: " << it << "\n";
          Rcout << "Acceptance Rates: " << acc_rates << "\n";
        }
        
        // reset within-cycle parameters
        it = 0;
        std::fill(xi_samp_counts.begin(),xi_samp_counts.end(),0);
        std::fill(acc_rates.begin(),acc_rates.end(),0.0);
      }
    }
    if(!quiet){
      Rcout << "\n...done!\n";
    }

    /* Main Sampling Phase */
    
    // allocate storage
    arma::mat chain_Xi(P*K,nstore);
    arma::mat chain_Xi0(P,nstore);
    arma::mat chain_II(K,nstore);
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
      for(int k = 0; k != K; ++k){
        Xik = Xi.col(k);
        Xik_prop = Xik + arma::trans(R.slice(k)) * arma::randn(P);
        double log_acc_prob = log_posterior_single(Xik_prop,Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),II(k),N(k),psi) - 
          log_posterior_single(Xik,Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),II(k),N(k),psi);
        if(R::runif(0,1) < std::exp(log_acc_prob)){
          Xi.col(k) = Xik_prop;
        }
      }
      
      // global params with gibbs proposals
      if(!sr_style){
        update_xi0(Xi,Xi0,V,V0); // modifies Xi0 in place
      }
      
      // update II
      update_ii(Xi,Y,Design_Matrices,gamma,T_1,II,N,psi,
                ii_prior_rate,ii_proposal_sd,1.0);
      
      update_Vparam(Xi,Xi0,Vparam,IGSR,lambda,ncovs,nbases,vparam_key);
      expand_Vparam(V,Vparam,lambda,ncovs,nbases);
      
      if(!(it < warmup) && !(it % thin)){
        int s = (it - warmup)/thin;
        chain_Xi.col(s) = Xi.as_col();
        chain_Xi0.col(s) = Xi0;
        chain_II.col(s) = II;
        chain_V.col(s) = Vparam;
      }
    }
    Rcout << " 100\n...done!\n\n";
    List chain_output;
    chain_output["Xi"] = chain_Xi;
    chain_output["Xi0"] = chain_Xi0;
    chain_output["II"] = chain_II;
    chain_output["V"] = chain_V;
    MCMC_output[chn] = chain_output;
  }
  return(MCMC_output);
}
