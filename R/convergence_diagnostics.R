convergence_diagnostics <- function(samples_matrix){
  if(!is.matrix(samples_matrix)){
    stop("Piss off!")
  }
  if(nrow(samples_matrix) %% 2){ # drop final row if total is odd
    samples_matrix <- samples_matrix[1:(nrow(samples_matrix) - 1),]
  }
  twoN <- nrow(samples_matrix)
  N <- twoN/2
  samples_matrix <- cbind(samples_matrix[1:N,],samples_matrix[(1+N):twoN,])
  P <- ncol(samples_matrix)
  
  grand_mean <- mean(c(samples_matrix))
  chain_means <- colMeans(samples_matrix)
  chain_var <- apply(samples_matrix,2,var)
  
  B <- N * sum((chain_means - grand_mean)^2) / (P - 1) # between chain variance
  W <- mean(chain_var) # within chain variance
  
  A <- (N - 1)*W/N + B/N; # aggregated variance estimate
  Rhat <- sqrt(A/W)
  
  Vt <- matrix(0,P,N)
  for(p in 1:P){
    for(n in 1:(N - 1)){
      Vt[p,n] = mean((samples_matrix[(n + 1):N,p] - samples_matrix[1:(N - n),p])^2)
    }
  }
  Vt <- colMeans(Vt)
  ac <- 1 - Vt/(2*A)
  st <- N - (N %% 2);
  BigP <- ac[seq(1,st - 1,by = 2)] + ac[seq(2,st,by = 2)]
  if(any(BigP < 0)){
    negBigPs <- which(BigP < 0)
    maxnn <- negBigPs[1] - 1 # maximum nonnegative
  }else{
    maxnn <- length(BigP) - 1
  }
  tau <- 1 + 2*sum(ac[1:(2*maxnn)])
  Neff <- round(P*N/tau)
  
  return(c(Rhat = Rhat, ESS = Neff))
}
