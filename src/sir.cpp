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

double log_negb(const NumericVector y, 
                  const NumericVector mu, 
                  const double DISP)
{
  double sz = 1/DISP;
  
  double res = 0;
  for(int j = 0; j < y.length(); ++j){
    res += R::dnbinom_mu(y[j],sz,mu[j],1);
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
                               const double IIP, // region's initial infectious population
                               const int N, // region population
                               const NumericVector vaccinations
){
  int J = beta.size();
  NumericVector nu(J); // initialize to length J vector of 0s
  
  // first case detection occurs in the interval t_0 -> t_1
  int offset = T_1 - 1; // translate index to Cpp: nu_1 is labeled nu_0 here.
  nu(offset) = IIP/N; // proportion of population infected at outbreak
  double state[] = {1 - nu(offset), nu(offset)}; // state at T_1
  double infections = 0, removals = 0, unvaccinated = 1;
  for(int j = offset + 1; j != J; ++j){
    infections = max(min(beta(j)*state[0]*state[1],
                         state[0] - (state[0]/unvaccinated)*vaccinations[j]),0);
    removals = max(min(gamma*state[1], state[1]),0);
    state[0] = state[0] - (state[0]/unvaccinated)*vaccinations[j] - infections;
    state[1] = state[1] + infections - removals;
    unvaccinated = max(unvaccinated - vaccinations[j],state[0]); // state[0] is a subset of unvaccinated and is nonnegative

    nu[j] = infections;
  }
  nu[nu == 0] = 1e-16;
  return N*(nu + abs(nu))/2; // get positive part to remove round-off negativity 
}

//[[Rcpp::export]]
NumericVector solve_susceptible(const NumericVector beta, // time-varying transmission rate
                               const double gamma, // inverse of average infectious period
                               const int T_1,   // time index of region's outbreak
                               const double IIP, // region's initial infectious population
                               const int N, // region population
                               const NumericVector vaccinations
){
  int J = beta.size();
  NumericVector susceptible(J, 1.0); // initialize to length J vector of 0s
  
  // first case detection occurs in the interval t_0 -> t_1
  int offset = T_1 - 1; // translate index to Cpp: nu_1 is labeled nu_0 here.
  double state[] = {1 - IIP/N, IIP/N}; // state at T_1
  double infections = 0, removals = 0, unvaccinated = 1;
  susceptible[offset] = state[0];
  for(int j = offset + 1; j != J; ++j){
    infections = max(min(beta(j)*state[0]*state[1],
                         state[0] - (state[0]/unvaccinated)*vaccinations[j]),0);
    removals = max(min(gamma*state[1], state[1]),0);
    state[0] = state[0] - (state[0]/unvaccinated)*vaccinations[j] - infections;
    state[1] = state[1] + infections - removals;
    unvaccinated = max(unvaccinated - vaccinations[j],state[0]); // state[0] is a subset of unvaccinated and is nonnegative
    susceptible[j] = state[0];
  }
  return susceptible; // get positive part to remove round-off negativity 
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
    // if(j < 70){
    //   for(int m = 0; m != u; ++m){
    //     delta[j] += nu[j - m] * psi[m];
    //   }
    // }else{
    //   for(int m = 0; m != u; ++m){
    //     delta[j] += nu[j - m] * psi[m]/(3.6);
    //   }
    // }
  }
  return delta;
}


// single region's contribution to the log posterior density
double log_posterior_single(const arma::vec & Xi, // P vector of region's parameters
                     const arma::vec & Xi_adjusted, // P vector adjusted to unrestricted values
                     const arma::vec & Xi0,     // P vector of global parameters
                     const arma::vec & V,       // variance parameters
                     const NumericVector Y, // region's event time series
                     const arma::mat dmat, // region's design matrices
                     const double gamma, // inverse average infectious period
                     const int T_1, // region's outbreak index
                     const double IIP, // region's initial infectious population
                     const double expected_iip, // prior mean IIP
                     const double DISP, // 'prior sample size' (1/nb dispersion)
                     const double expected_disp, // prior mean DISP
                     const int N, // region population
                     const NumericVector psi, // interval event probability
                     const NumericVector vaccinations,
                     const int discount_period_length, // length of the "discount" period at the beginning of the outbreak
                     const double discount_period_disp
){

  arma::vec B = arma::exp( dmat * Xi );
  // // check if parameters are within prior support first
  // bool B_out_of_bounds = arma::any(arma::vectorise(B) < 0);
  // if(B_out_of_bounds || IIP > N){ // B or IIP out of bounds?
  //   return(std::log(0)); // prior(beta) = 0
  // }
  
  double res = 0;
  // add log prior contribution (Xi_adjusted)
  arma::mat X(Xi_adjusted - Xi0);
  res += -0.5 * arma::as_scalar( arma::dot(((X*(X.t())).eval()).diag(), 1/V) );
  
  // int P = Xi.n_rows;// temp 
  // arma::vec mu0(P,arma::fill::zeros); // temp
  // mu0(0) = 2; // temp
  // arma::vec post_mean = post_var * ((K/V) % mean(Xi,1) + (1/V0) % mu0); // temp edited // rowMeans(Xi)
  
  
  // add log prior contribution (IIP)
  // prior on IIP is exponential with mean `expected_iip`
  res += -IIP/expected_iip;

  // prior on DISP is Normal+ (Normal folded at 0) with mean `expected_disp`
  res += -(DISP*DISP/(expected_disp*expected_disp*3.1415926536));   
  
  // add log likelihood contribution
  NumericVector events = solve_events(solve_infections(as<NumericVector>(wrap(B)),gamma,T_1,IIP,N,vaccinations),psi);

  if(discount_period_length == 0){
    res += log_negb(Y,events,DISP);
  }else{
    // subsetting machinery
    size_t J = Y.length();
    NumericVector sel1(discount_period_length);
    NumericVector sel2(J - discount_period_length);
    for(int j = 0; j < J; ++j){
      if(j < discount_period_length){
        sel1[j] = j;
      }else{
        sel2[j - discount_period_length] = j;
      }
    }
  
    NumericVector Yhead = Y[sel1];
    NumericVector Ytail = Y[sel2];
    NumericVector Ehead = events[sel1];
    NumericVector Etail = events[sel2];
  
    res += log_negb(Yhead, Ehead, discount_period_disp); 
    res += log_negb(Ytail, Etail, DISP);
  }
  return(res);
}

void update_xi0(arma::mat const & Xi_adjusted,
                arma::vec & Xi0,
                arma::vec const & V,
                arma::vec const & V0
)
{
  int K = Xi_adjusted.n_cols;
  // int P = Xi.n_rows;// temp 
  // arma::vec mu0(P,arma::fill::zeros); // temp
  // mu0(0) = 2; // temp
  arma::mat post_var = arma::diagmat(1/(K*(1/V) + (1/V0))); // element-wise inverse
  // arma::vec post_mean = post_var * ((K/V) % mean(Xi,1) + (1/V0) % mu0); // temp edited // rowMeans(Xi)
  arma::vec post_mean = post_var * (K/V) % mean(Xi_adjusted,1); // rowMeans(Xi)
  Xi0 = mvrnorm(post_mean, post_var);
}

void update_Vparam(arma::mat const & Xi_adjusted,
                   arma::vec const & Xi0,
                   arma::vec & Vparam,
                   NumericMatrix IGSR, // Inverse-Gamma Shape and Rate hyperparams
                   arma::vec const & lambda,
                   int const ncovs,
                   int const nbases,
                   bool const vparam_key[]
)
{
  arma::mat X(Xi_adjusted.each_col() - Xi0);
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
                 const List TildeOff, // List of K matrices (removes tilde from coefficient set k)
                 const NumericMatrix vaccinations, // JxK matrix of vaccination data
                 const arma::vec & lambda, // GP basis eigenvalues
                 const arma::vec & V0param, // pre-expanded variance hyperparameters
                 const NumericMatrix IGSR, // (3x2 if hierarchical, 1x2 ow) Gamma Shape and Rate hyperparameters
                 const double gamma, // inverse infectious period
                 const IntegerVector T_1, // outbreak indices
                 const double expected_iip, // expected infectious population
                 const double expected_disp,
                 const IntegerVector N, // region populations
                 const NumericMatrix psi, // MxK interval event probability
                 const int discount_period_length,
                 const double discount_period_disp,
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
    arma::mat Xi_adjusted(P,K);
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
      arma::mat tildeoff_k = TildeOff[k];
      Xi_adjusted.col(k) = tildeoff_k * Xik;
      if(!quiet){
        Rcout << " " << k + 1;
      }
    }
    arma::vec Xik_prop = Xik; // armadillo does a deep copy
    arma::vec Xik_prop_adjusted = Xik; // these values are just placeholders

    // initialize initial infectious population
    arma::rowvec IIP(K);
    for(int k = 0; k != K; ++k){
      IIP[k] = R::rgamma(1,expected_iip); // randomly initialize IIP ~ exp(1/expected_iip)
    }
    double log_IIPk = 0; // initial value doesn't matter
    double log_IIPk_prop = 0; // initial value doesn't matter
    
    // initialize 'prior sample size' (inverse of nb dispersion parameter)
    arma::rowvec DISP(K);
    for(int k = 0; k != K; ++k){
      DISP[k] = R::rgamma(1,1) + 1;
    }
    double log_DISPk = 0; // initial value doesn't matter
    double log_DISPk_prop = 0; // initial value doesn't matter
    
    // Join IIP and then DISP to Xi as new rows
    arma::mat Xi_and_log_IIP = arma::join_cols(Xi,arma::log(IIP));
    Xi_and_log_IIP = arma::join_cols(Xi_and_log_IIP,arma::log(DISP));
    arma::vec Xik_and_log_IIPk_prop(P+2);
    
    
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
    //Xi0(0) = 3.5;
    if(!sr_style){
      update_xi0(Xi,Xi0,V,V0); // go ahead and initialize Xi0 among the Xi
    }

    if(!quiet){
      Rcout << " done!\n";
      Rcout << "Begin adaptation... \n\t";
    }
    
    /* Adaptation Phase */
    
    // prep to adapt Xi proposal covariance
    arma::vec C0diag(P + 2,arma::fill::ones);
    C0diag.subvec(1 + ncovs, P - 1) = lambda/10;
    arma::mat C0 = arma::diagmat(C0diag);
    arma::cube R(P + 2, P + 2, K);
    for(int k = 0; k != K; ++k){
      R.slice(k) = arma::chol(C0);
   }
    
    // arguments ncycles, samps per cycle, quiet
    IntegerVector xi_samp_counts(K,0);
    IntegerVector xi_samps_since_last_accepted(K,0);
    NumericVector acc_rates(K,0.0);
    int completed_cycle_count = 0;
    bool adapting = true;
    int it = 0; // within-cycle iteration count
    arma::cube ap_Xi_and_log_IIP(samps_per_cycle,P + 2,K);
    while(adapting){
      it += 1; // update iteration number
      if((it % 1000) == 0){checkUserInterrupt();}
      for(int k = 0; k != K; ++k){
        if(xi_samp_counts[k] < samps_per_cycle){
          Xik_and_log_IIPk_prop = Xi_and_log_IIP.col(k) + arma::trans(R.slice(k)) * arma::randn(P + 2);
          Xik_prop = Xik_and_log_IIPk_prop.subvec(0,P - 1);
          arma::mat tildeoff_k = TildeOff[k];
          Xik_prop_adjusted = tildeoff_k * Xik_prop;
          log_IIPk_prop = Xik_and_log_IIPk_prop[P]; log_DISPk_prop = Xik_and_log_IIPk_prop[P + 1];
          Xik = Xi_and_log_IIP(arma::span(0,P - 1),k); log_IIPk = Xi_and_log_IIP(P,k); log_DISPk = Xi_and_log_IIP(P+1,k);
          double log_acc_prob = log_posterior_single(Xik_prop,Xik_prop_adjusted,Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),std::exp(log_IIPk_prop),expected_iip,std::exp(log_DISPk_prop),expected_disp,N(k),psi(_,k),vaccinations(_,k),discount_period_length,discount_period_disp) - 
            log_posterior_single(Xik,Xi_adjusted.col(k),Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),std::exp(log_IIPk),expected_iip,std::exp(log_DISPk),expected_disp,N(k),psi(_,k),vaccinations(_,k), discount_period_length, discount_period_disp) + 
            (log_IIPk_prop - log_IIPk) + (log_DISPk_prop - log_DISPk); // add Jacobians due to variable transforms
          if(R::runif(0,1) < std::exp(log_acc_prob)){
            Xi_and_log_IIP.col(k) = Xik_and_log_IIPk_prop;
            ap_Xi_and_log_IIP(xi_samp_counts[k],0,k,arma::size(1,P+2,1)) = Xik_and_log_IIPk_prop;
            xi_samp_counts[k] += 1;
            xi_samps_since_last_accepted[k] = 0;
            acc_rates[k] = (1 + (it - 1)*acc_rates[k])/it;
          }else{
            if(completed_cycle_count > ncycles){ // stop being greedy
              ap_Xi_and_log_IIP(xi_samp_counts[k],0,k,arma::size(1,P+2,1)) = Xi_and_log_IIP.col(k);
              xi_samp_counts[k] += 1;
            }
            xi_samps_since_last_accepted[k] += 1;
            acc_rates[k] = (0 + (it - 1)*acc_rates[k])/it;
          }
          if(xi_samps_since_last_accepted[k] > 100){ //completed_cycle_count < ncycles && 
            R.slice(k) = R.slice(k)/3.16; // divide proposal covariance by 10
            xi_samps_since_last_accepted[k] = 0;
          }
        }
      }
      Xi = Xi_and_log_IIP.rows(0,P - 1);
      for(int k = 0; k < K; ++k){
        arma::mat tildeoff_k = TildeOff[k];
        Xi_adjusted.col(k) = tildeoff_k * Xi.col(k);
      }
      IIP = arma::exp(Xi_and_log_IIP.row(P));
      DISP = arma::exp(Xi_and_log_IIP.row(P+1));
      //Rcout << "Finished XI update \n";
      // global params with gibbs proposals
      if(!sr_style){
        update_xi0(Xi_adjusted,Xi0,V,V0); // modifies Xi0 in place
      }
      //Rcout << "Finished XI0 update \n";
      // update variance parameters
      update_Vparam(Xi_adjusted,Xi0,Vparam,IGSR,lambda,ncovs,nbases,vparam_key);
      expand_Vparam(V,Vparam,lambda,ncovs,nbases);
      //Rcout << "Finished  V update \n";
      // is the adaptation cycle complete?
      if(is_true(all(xi_samp_counts == samps_per_cycle))){
        // increment the completed cycle count
        completed_cycle_count += 1;
        
        // do we need to keep adapting? (are we short on cycles, or are acc_rates unacceptable)
        adapting = ((completed_cycle_count < ncycles) || is_true(any(acc_rates < 0.3)) || is_true(any(acc_rates > 0.4))) && (completed_cycle_count < ncycles + 10);
        
        // compute new Xi proposal covariances
        for(int k = 0; k!= K; ++k){
          // magic number comes from 2.38^2, see article Haario et al (2001)
          if(completed_cycle_count < ncycles){
            //Rcout << k << '\n';
	    arma::mat SigHat = 5.66*arma::cov(ap_Xi_and_log_IIP.slice(k))/(P + 2);
	    // Take the robust cholesky decomposition
	    bool chol_success = false;
	    while(chol_success == false){
	      chol_success = arma::chol(R.slice(k), SigHat);
	      if(chol_success == false){
		SigHat += C0 * 1e-6;
	      }
	    }
          }else{ // after ncycle, begin retaining covariance information across cycles
            int extra_cycle_count = completed_cycle_count - (ncycles - 1);
            arma::mat oldSigHat = arma::trans(R.slice(k)) * R.slice(k);
            arma::mat newSigHat = 5.66*arma::cov(ap_Xi_and_log_IIP.slice(k))/(P + 2);
            arma::mat meanSigHat = ((extra_cycle_count*samps_per_cycle - (P + 2))*oldSigHat + 
              (samps_per_cycle - (P + 2))*newSigHat)/((extra_cycle_count+1)*samps_per_cycle - (P + 2));
	    // take the robust cholesky decomposition
	    bool chol_success =	false;
            while(chol_success == false){
              chol_success = arma::chol(R.slice(k), meanSigHat);
	      if(chol_success == false){
                meanSigHat += C0 * 1e-6;
              }
	    }
          }
        }

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
    arma::mat chain_IIP(K,nstore);
    arma::mat chain_DISP(K,nstore);
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
        Xik_and_log_IIPk_prop = Xi_and_log_IIP.col(k) + arma::trans(R.slice(k)) * arma::randn(P + 2);
        Xik_prop = Xik_and_log_IIPk_prop.subvec(0,P - 1);
        arma::mat tildeoff_k = TildeOff[k];
        Xik_prop_adjusted = tildeoff_k * Xik_prop;
        log_IIPk_prop = Xik_and_log_IIPk_prop[P]; log_DISPk_prop = Xik_and_log_IIPk_prop[P + 1];
        Xik = Xi_and_log_IIP(arma::span(0,P - 1),k); log_IIPk = Xi_and_log_IIP(P,k); log_DISPk = Xi_and_log_IIP(P+1,k);
        double log_acc_prob = log_posterior_single(Xik_prop,Xik_prop_adjusted,Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),std::exp(log_IIPk_prop),expected_iip,std::exp(log_DISPk_prop),expected_disp,N(k),psi(_,k),vaccinations(_,k),discount_period_length,discount_period_disp) - 
          log_posterior_single(Xik,Xi_adjusted.col(k),Xi0,V,Y(_,k),Design_Matrices[k],gamma,T_1(k),std::exp(log_IIPk),expected_iip,std::exp(log_DISPk),expected_disp,N(k),psi(_,k),vaccinations(_,k),discount_period_length,discount_period_disp) + 
          (log_IIPk_prop - log_IIPk) + (log_DISPk_prop - log_DISPk); // add Jacobians due to variable transforms;
        if(R::runif(0,1) < std::exp(log_acc_prob)){
          Xi_and_log_IIP.col(k) = Xik_and_log_IIPk_prop;
        }
      }
      Xi = Xi_and_log_IIP.rows(0,P - 1);
      for(int k = 0; k < K; ++k){
        arma::mat tildeoff_k = TildeOff[k];
        Xi_adjusted.col(k) = tildeoff_k * Xi.col(k);
      }
      IIP = arma::exp(Xi_and_log_IIP.row(P));
      DISP = arma::exp(Xi_and_log_IIP.row(P+1));
      
      // global params with gibbs proposals
      if(!sr_style){
        update_xi0(Xi_adjusted,Xi0,V,V0); // modifies Xi0 in place
      }
      
      update_Vparam(Xi_adjusted,Xi0,Vparam,IGSR,lambda,ncovs,nbases,vparam_key);
      expand_Vparam(V,Vparam,lambda,ncovs,nbases);
      
      if(!(it < warmup) && !(it % thin)){
        int s = (it - warmup)/thin;
        chain_Xi.col(s) = Xi.as_col();
        chain_Xi0.col(s) = Xi0;
        chain_IIP.col(s) = IIP.as_col();
        chain_DISP.col(s) = DISP.as_col();
        chain_V.col(s) = Vparam;
      }
    }
    Rcout << " 100\n...done!\n\n";
    List chain_output;
    chain_output["Xi"] = chain_Xi;
    chain_output["Xi0"] = chain_Xi0;
    chain_output["IIP"] = chain_IIP;
    chain_output["DISP"] = chain_DISP;
    chain_output["V"] = chain_V;
    MCMC_output[chn] = chain_output;
  }
  return(MCMC_output);
}
