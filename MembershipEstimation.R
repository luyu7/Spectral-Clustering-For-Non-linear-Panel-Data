library(rlist)
library(quantreg)
library(base)
library(expm)
library(dplyr)
library(combinat)
library(stats)
library(gplots)
library(wordspace)
library(cluster)
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
  
  sapply(sub, function(i) as.numeric(seq((i-1)*L+1,i*L,1)))
}

variation.inf <- function( gammaA, gammaB, covA, covB){
  # Calclulate the dissimilairty between subject A and subject B. 
  #
  # Args:
  #      gammaA: coefficient estimates for subject A
  #      gammaB: coefficient estimates for subject B
  #      covA: variance estimator of gammaA
  #      covB: variance estimator of gammaB
  #
  # Returns:
  #      Dissimlairty between A and B
  
  inv.cov <- solve(covA + covB)
  scaled.variation <- ( sqrtm(inv.cov) %*% (gammaA-gammaB))
  result <- norm(scaled.variation, type = "I")
  return(result)
}

membership.est <- function(X, y, n, L, G)
{
  # Estimate the membership using the proposed spectral clustering method.
  #
  # Args:
  #      X: (n*T) * p design matrix
  #      y: (n*T * 1) response
  #      n: number of subjects
  #      L: number of observations for each subject
  #      G: number of groups
  #
  # Returns:
  #      Estimated membership
  
  membership0 <- as.list(seq(1,n))
  gamma.list <- list()
  beta.list <- list()
  cov.list <- list()
  sub.cov.list <- list()
  for(k in 1:n)
  {
    ind.temp <- as.vector(reconstruct(membership0[[k]], L))
    new.X <- as.matrix(X[ind.temp, ])
    new.y <- as.matrix(y[ind.temp])
    mod <- rq(new.y ~ new.X , tau = tau, method = "fn")
    gamma.list[[k]] <- as.numeric(mod$coefficients)
    beta.list[[k]] <- gamma.list[[k]][2:(p+1)]
    cov.list[[k]] <- summary.rq(mod, se = "nid", covariance = T)$cov
    sub.cov.list[[k]] <- cov.list[[k]][-1,-1]
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
  
  # Construct the diagonal matrix.
  D <- diag(apply(A, 1, sum)) 
  
  # Compute a normalized graph Laplacian.
  Lap <- solve(sqrtm(D)) %*% (D-A)  %*% solve(sqrtm(D))
  
  # Find the k eigenvectors of Lap associated to the k smallest eigenvalues, 
  # and form the matrix Z by stacking the k eigenvectors in columns.
  eigenv <- eigen(Lap, symmetric=TRUE)
  Z <- eigenv$vectors[,(ncol(eigenv$vectors)-G+1):ncol(eigenv$vectors)]
  
  # Normalize rows of Z.
  Z <- normalize.rows(Z)
  
  # Standard k-means clustering. 
  km <- kmeans(Z, centers=G, nstart = 20, algorithm = "Hartigan-Wong", iter.max = 30)
  result <- km$cluster
  return(result)
}
