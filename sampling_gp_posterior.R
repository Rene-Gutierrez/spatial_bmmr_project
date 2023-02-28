################################################################################
# Draws Samples from a Gaussian Process
################################################################################

################################################################################
# Settings
# Number of Samples
sam <- 3
# Samples the Observations
source("./gaussian_process_simulation.R")
################################################################################

################################################################################
# New Regular Grid
# Number of Points
NG <- 100
# Grid
SG <- matrix(data = NA, nrow = NG, ncol = D)
for(d in 1:D){
  SG[, d] <- seq(from = 0, to = grid_size[d], length.out = NG)
}
################################################################################

################################################################################
# Distance Matrices
# Distance Matrix for the New Grid
DMnn <- matrix(data = NA, nrow = NG, ncol = NG)
for(n in 1:NG){
  for(nn in 1:NG){
    DMnn[n, nn] <- sum((SG[n, ] - SG[nn, ])^2)
  }
}
# Distance Matrix Between the the New Grid and the Observations
DMno <- matrix(data = NA, nrow = NG, ncol = N)
for(n in 1:NG){
  for(nn in 1:N){
    DMno[n, nn] <- sum((SG[n, ] - S[nn, ])^2)
  }
}
################################################################################

################################################################################
# Samples f
# Auxilary Variables
Sno <- t^2 * exp(- DMno / (2 * l^2))
Snn <- t^2 * exp(- DMnn / (2 * l^2))
# Samples f
mea <- Sno %*% solve(Soo + s2 * diag(nrow = N), y)
Sig <- Snn - Sno %*% solve(Soo + s2 * diag(nrow = N), t(Sno))
fs  <- mvtnorm::rmvnorm(n = sam, mean = mea, sigma = Sig)
# Samples s2
################################################################################


################################################################################
# Plots
# Plotting Dimensions
ymax <-   2
ymin <- - 2
# Plots the Observations
plot(x    = S,
     y    = y,
     xlim = c(0, grid_size),
     ylim = c(ymin, ymax),
     pch  = 18)
# Plots the Gaussian Process Draws
for(i in 1: sam){
  par(new = TRUE)
  plot(x    = SG,
       y    = fs[i, ],
       xlim = c(0, grid_size),
       ylim = c(ymin, ymax),
       type = 'l',
       col  = rgb(0, 0, 1, 0.25))
}
par(new = TRUE)
plot(x    = S,
     y    = y,
     xlim = c(0, grid_size),
     ylim = c(ymin, ymax),
     pch  = 18)
################################################################################