library(MASS)
library(rlist)
library(quantreg)
library(base)
library(expm)
library(dplyr)
library(combinat)
library(stats)
options(digits = 8)

rhotau <- function(x, tau){
  # Check function.
  #
  # Args:
  #      x: input
  #      tau: quantile
  #
  # Returns:
  #      The value of check function with input x and quanitle tau
  
  (tau - 1*(x <= 0)) * x
}


reconstruct <- function(sub, L){
  # Generate the indices for a subject.
  #
  # Args:
  #      sub: index of the subject
  #        L: number of observations for the subject
  #
  # Returns:
  #      Indices of the observations for the given subject 
  
  sapply(sub, function(i) as.numeric(seq((i-1)*L+1, i*L, 1)))
}

variation.inf <- function( gammaA, gammaB, covA, covB){
  # Calclulate the scaled dissimilarity between subject A and subject B. 
  #
  # Args:
  #      gammaA: coefficient estimates for subject A
  #      gammaB: coefficient estimates for subject B
  #      covA: variance estimator of gammaA
  #      covB: variance estimator of gammaB
  #
  # Returns:
  #      Scaled dissimilarity between A and B
  
  scaled.variation <- expm((-0.5) * logm(covA + covB)) %*% (gammaA - gammaB)
  result <- norm(scaled.variation, type = "I")
  result <- result * 2 / sqrt(log(L) * log(n))
  return(result)
}

G.est <- function(X, y, n, L, Gmax){
  # Estimate the number of groups.
  #
  # Args: 
  #      X: (n*T) * p design matrix
  #      y: (n*T * 1) response
  #      n: number of subjects
  #      L: number of observations for each subject
  #      Gmax: maximum number of clusters to consider
  #
  # Returns:
  #      Estimated number of groups
  
  membership0 <- as.list(seq(1, n))
  gamma.list <- list()
  beta.list <- list()
  cov.list <- list()
  sub.cov.list <- list()
  for(k in 1:n)
  {
    set.seed(1)
    ind0 <- as.vector(reconstruct(membership0[[k]], L))
    new.X <- as.matrix(X[ind0, ])
    new.y <- as.matrix(y[ind0])
    mod <- rq(new.y ~ new.X, tau = tau, method = "fn")
    beta.list[[k]] <- as.numeric(mod$coefficients)[2:3]
    sub.cov.list[[k]] <- summary.rq(mod, se = "nid", covariance = T)$cov[-1,-1]
  }
  
  # Construct the empirical dissimilarity matrix.
  variation.mat <- matrix(NA, n, n)
  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    {
      variation.mat[i,j] <- variation.inf(beta.list[[i]], beta.list[[j]], sub.cov.list[[i]], sub.cov.list[[j]])
    }
  }
  
  # Construct the empirical adjacency matrix.
  A <- variation.mat
  diag(A) <- 0
  A[lower.tri(A)] = t(A)[lower.tri(A)]
  A <- exp(-A)
  dv <- 1/sqrt(rowSums(A))
  l <- dv * A %*% diag(dv)
  decomp <- eigen(l)
  evals <- as.numeric(decomp$values)
  
  # Estimate the number of groups using the proposed statistics (see eq (2) in the paper).
  diffs <- diff(evals)
  est.G <- which.max(abs(diffs)[1:maxG]/(evals[2:(maxG+1)]))
  return(est.G)
}