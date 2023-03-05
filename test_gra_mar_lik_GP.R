################################################################################
# Test Gradient Marginal Likelihood Gaussian Process
################################################################################

# Library
source("./gaussian_process_simulation.R")
source("./gra_mar_lik_GP.R")

# Auxiliary Variables
K <- t2 * exp(- DM / (2 * l2))
C <- K + s2 * diag(nrow = N)
h <- 0.00000001

# Runs the Gradient
out1 <- gra_mar_lik_GP(s2 = s2,
                       t2 = t2,
                       l2 = l2,
                       y  = y,
                       DM = DM)

# Test for s2
Kh   <- t2 * exp(- DM / (2 * l2))
Ch   <- K + (s2 + h) * diag(nrow = N)
out2 <- ((-log(det(Ch)) / 2 - t(y) %*% solve(Ch, y) / 2) - (-log(det(C)) / 2 - t(y) %*% solve(C, y) / 2)) / h

# Test for t2
Kh   <- (t2 + h) * exp(- DM / (2 * l2))
Ch   <- Kh + s2 * diag(nrow = N)
out3 <- ((-log(det(Ch)) / 2 - t(y) %*% solve(Ch, y) / 2) - (-log(det(C)) / 2 - t(y) %*% solve(C, y) / 2)) / h

# Test for t2
Kh   <- t2 * exp(- DM / (2 * (l2 + h)))
Ch   <- Kh + s2 * diag(nrow = N)
out4 <- ((-log(det(Ch)) / 2 - t(y) %*% solve(Ch, y) / 2) - (-log(det(C)) / 2 - t(y) %*% solve(C, y) / 2)) / h
