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
N <- 200
# Variance
s2  <- 1 / 2
t2  <- 2
l2  <- 1 / 2
################################################################################

################################################################################
# Point Locations
S <- matrix(data = NA, nrow = N, ncol = D)
for(n in 1:N){
  for(d in 1:D){
    S[n, d] <- runif(n = 1, max = grid_size[d])
  }
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
################################################################################

################################################################################
# Simulates
# Covariance Matrix
Soo <- t2 * exp(- DM / (2 * l2))
# Simulates Observations
y   <- as.vector(mvtnorm::rmvnorm(n = 1, sigma = Soo + s2 * diag(nrow = N)))
################################################################################

################################################################################
# Plots
if(D == 1){
  # Plotting Dimensions
  ymax <- max(   2)
  ymin <- min(  -2)
  # Plotting Order
  o    <- order(S[,1])
  # Plots
  plot(x    = S[,1][o],
       y    = y[o],
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
