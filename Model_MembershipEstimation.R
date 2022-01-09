rm(list=ls())
source("MembershipEstimation.R")



# Data generation via Model 2 with n=30, T=15, tau=0.5, and t(3) error.
n <- 30
L <- 15
tau <- 0.5
G <- 3
d <- matrix(rep(1:G, rep(1/G,G)*n), n, 1)
truegroups <- lapply(permn(1:G), 
                     function(x) as.numeric(levels(factor(d,levels = c(1,2,3),labels=x)))[factor(d,levels = c(1,2,3),labels=x)])
ind <- rep(d, each=L)
bmatrix <- matrix(c(-5,1,0,1,3,1),nrow = 2,ncol = 3)
alpha <- rep(1,n)
intercept <- rep(alpha, each=L)
set.seed(1)
gamma.all <- rnorm(3*n)
gamma1 <- rep(gamma.all[1:n], rep(L,n))
gamma2 <- rep(gamma.all[(n+1):(2*n)], rep(L,n))
set.seed(1)
x1 <- 0.5*intercept + gamma1 + rnorm(n*L,0,sd=1)
x2 <- 0.5*intercept + gamma2 + rnorm(n*L,0,sd=sqrt(0.05))
X <- cbind(x1,x2)
p <- ncol(X)

# Number of iterations for the simulation.
R <- 100

# Store the estimated membership.
record <- list()  

for (r in 1:R)
{
  set.seed(r)
  e <- rt(n*L,df=3)
  e <- as.vector(t(e))
  y <- rep(0,n*L)
  for(g in 1:G)
  {
    y[ind==g] <- intercept[ind==g] + X[ind==g,] %*% bmatrix[,g] + e[ind==g]
  }
  record[[r]] <- membership.est(X, y, n, L, G)
}

# Compute perfect and average match percentage.
count.perfect <- c()
count.ave <- c()
for(i in 1:R)
{
  temp.membership <- record[[i]]
  count.perfect[i] <- any(sapply(truegroups, function(x, want) isTRUE(all.equal(x, want)), temp.membership))
  count.ave[i] <- max(sapply(truegroups, function(x, want) sum(want==x)/n, temp.membership))
}

# Print the results.
cat("---  Spectral Clustering  ---", "Perfect Match:", sum(count.perfect)/R, " Average Match:", mean(count.ave))

