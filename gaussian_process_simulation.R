################################################################################
# Simulates A Gaussian Process 
################################################################################

################################################################################
# Settings
# Libraries

# Dimensions
D <- 1
# Grid Size
grid_size <- rep(10, D)
# Number of Points
K <- 25
M <- 100
N <- 1000
# Variance
s2  <- 1^2
t2  <- (1 / 3)^2
l2  <- (1 / 4)^2
################################################################################

################################################################################
# Point Locations
S <- matrix(data = NA, nrow = N, ncol = D)
for(n in 1:N){
  for(d in 1:D){
    S[n, d] <- runif(n = 1, max = grid_size[d])
  }
}
Sk <- matrix(data = NA, nrow = M, ncol = D)
for(d in 1:D){
  Sk[, d] <- seq(0, grid_size[d], length.out = M)
}
################################################################################

################################################################################
# Distance Matrix
DM <- matrix(data = 0, nrow = N, ncol = N)
for(n in 1:N){
  for(nn in 1:N){
    DM[n, nn] <- sum((S[n, ] - S[nn, ])^2)
  }
}
DM  <- sqrt(DM)
# Dkk
Dkk <- matrix(data = 0, nrow = M, ncol = M)
for(m in 1:M){
  for(mm in 1:M){
    Dkk[m, mm] <- sum((Sk[m, ] - Sk[mm, ])^2)
  }
}
Dkk <- sqrt(Dkk)
# Dks
Dks <- matrix(data = 0, nrow = M, ncol = N)
for(m in 1:M){
  for(n in 1:N){
    Dks[m, n] <- sum((Sk[m, ] - S[n, ])^2)
  }
}
Dks <- sqrt(Dks)
################################################################################

################################################################################
# Simulates
# Covariance Matrix
Soo <- t2 * exp(- DM / l2)
# Simulates Observations
y <- matrix(data = NA, nrow = N, ncol = K)
w <- as.vector(mvtnorm::rmvnorm(n = 1, sigma = Soo))
for(k in 1:K){
  y[, k] <- w + rnorm(n = N, sd = sqrt(s2))
}
################################################################################

################################################################################
# Plots
if(D == 1){
  # Plotting Dimensions
  ymax <- max(   3)
  ymin <- min(  -3)
  # Plotting Order
  o    <- order(S[,1])
  # Plots
  plot(x    = S[,1][o],
       y    = w[o],
       ylim = c(ymin, ymax),
       type = 'l')
  for(k in 1:K){
    par(new=TRUE)
    plot(x    = S[,1][o],
         y    = y[o, k],
         ylim = c(ymin, ymax),
         col  = rgb(0, 0, 1, 0.1),
         type = 'l')
  }
  par(new=TRUE)
  plot(x    = S[,1][o],
       y    = w[o],
       ylim = c(ymin, ymax),
       type = 'l')
} else {
  # Color Limits
  cmin <-  (4 * s2)
  cmax <- -(4 * s2)
  ylev <- (y - cmin) / (cmax - cmin)
  # Plotting Dimensions
  xmax <- max(S[, 1])
  xmin <- min(S[, 1])
  ymax <- max(S[, 2])
  ymin <- min(S[, 2])
  # Plotting Area
  plot(NULL,
       type = 'l',
       ylim = c(ymin, ymax),
       xlim = c(xmin, xmax),
       xlab = "",
       ylab = "",
       # xaxt = 'n',
       # yaxt = 'n',
       main = "")
  # Plots
  rect(xleft   = S[, 1] - 0.01,
       ybottom = S[, 2] - 0.01,
       xright  = S[, 1] + 0.01,
       ytop    = S[, 2] + 0.01,
       col     = rgb(ylev, 0, 1 - ylev),
       border  = NA)
}
################################################################################
