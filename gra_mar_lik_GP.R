################################################################################
# Gradient Marginal Likelihood GP
################################################################################

gra_mar_lik_GP <- function(s2,
                           t2,
                           l2,
                           y,
                           DM){
  ##############################################################################
  # Dimensions
  N   <- length(y)
  ##############################################################################
  
  ##############################################################################
  # Gradient
  # Covariance Function
  K    <- t2 * exp(- DM / (2 * l2))
  eC   <- eigen(K + s2 * diag(nrow = N))
  iC   <- solve(K + s2 * diag(nrow = N))
  # Derivative s2
  g_s2 <-        sum(1 / eC$values) / 2
  g_s2 <- g_s2 - t(y) %*% eC$vectors %*% diag(1 / eC$values^2) %*% t(eC$vectors) %*% y / 2
  # Derivative t2
  g_t2 <-        sum(diag(solve(K + s2 * diag(nrow = N), K / t2))) / 2
  g_t2 <- g_t2 - t(y) %*% eC$vectors %*% diag(1 / eC$values^2) %*% t(eC$vectors) %*% K %*% y / (2 * t2)
  # Derivative t2
  g_l2 <-        sum(diag(solve(K + s2 * diag(nrow = N), (K * DM / (2 * l2^2)) ))) / 2
  g_l2 <- g_l2 - t(y) %*% iC %*% (K * DM / (2 * l2^2)) %*% iC %*% y / 2
  ##############################################################################
  
  ##############################################################################
  # Returns Values
  return(c(g_s2, g_t2, g_l2))
}