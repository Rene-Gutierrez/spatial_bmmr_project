################################################################################
# Gradient Transformed Marginal Likelihood GP
################################################################################

tra_gra_mar_lik_GP <- function(ls2,
                               lt2,
                               ll2,
                               y,
                               DM){
  ##############################################################################
  # Dimensions
  N   <- length(y)
  ##############################################################################
  
  ##############################################################################
  # Transforms the variables
  s2 <- exp(ls2)
  t2 <- exp(lt2)
  l2 <- exp(ll2)
  ##############################################################################
  
  ##############################################################################
  # Gradient
  # Covariance Function
  K    <- exp(lt2) * exp(- DM / (2 * exp(ll2)))
  C    <- K + exp(ls2) * diag(nrow = N)
  iC   <- solve(C)
  Cy   <- iC %*% y
  # Derivative s2
  g_ls2 <-       - exp(ls2) * sum(diag(iC)) / 2
  g_ls2 <- g_ls2 + t(Cy) %*% Cy * exp(ls2) / 2
  g_ls2 <- g_ls2 + 1
  # Derivative t2
  g_lt2 <-       - sum(diag(iC %*% K)) / 2
  g_lt2 <- g_lt2 + t(Cy) %*% K %*% Cy / 2
  g_lt2 <- g_lt2 + 1
  # Derivative t2
  g_ll2 <-       - sum(diag(iC %*% (K * DM))) / (2 * exp(ll2)) / 2
  g_ll2 <- g_ll2 + t(Cy) %*% (K * DM) %*% Cy  / (2 * exp(ll2)) / 2
  g_ll2 <- g_ll2 + 1
  ##############################################################################
  
  ##############################################################################
  # Returns Values
  return(c(g_ls2, g_lt2, g_ll2))
}