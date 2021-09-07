rm(list=ls())
source("NumGrpsEstimation.R")

# Data generation via Model 2 with n=30, T=15, tau=0.5, and t(3) error.
n <- 30
L <- 15
tau <- 0.5
G <- 3
set.seed(1)
d <- sample(c(1:G), size = n, replace = TRUE)
ind <- rep(d,each=L)
bmatrix <- matrix(c(-5,1,0,1,3,1), nrow = 2, ncol = 3)
set.seed(1)
alpha <- runif(n)
intercept <- rep(alpha, each=L)
set.seed(1)
gamma.all <- rnorm(3*n)
gamma1 <- rep(gamma.all[1:n], rep(L,n))
gamma2 <- rep(gamma.all[(n+1):(2*n)], rep(L,n))
set.seed(1)
x1 <- 0.5*intercept + gamma1 + rnorm(n*L, 0, sd=1)
x2 <- 0.5*intercept + gamma2 + rnorm(n*L, 0, sd=sqrt(0.05))
X <- cbind(x1,x2)
p <- ncol(X)

# Maximum number of clusters to consider.
maxG <- 10

# Number of iterations for the simulation.
R <- 100

# Store the percentage of estimated number of groups.
res.ratio <- rep(0,R)

for (r in 1:R){
  set.seed(r)
  e <- rt(n*L, df=3)
  e <- as.vector(t(e))
  y <- rep(0,n*L)
  for(g in 1:G)
  {
    y[ind==g] <- intercept[ind==g] + X[ind==g,] %*% bmatrix[,g] + e[ind==g]
  }
  res.ratio[r] <- G.est(X, y, n, L, Gmax)
} 

# Print the results.
prop.table( table(res.ratio) )

