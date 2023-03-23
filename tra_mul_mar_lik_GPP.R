################################################################################
# Marginal Likelihood Gaussian Predictive Process (Multiple and Transformed)
################################################################################

tra_mul_mar_lik_GPP <- function(s2,
                                ga,
                                ph,
                                yb,
                                Dkk,
                                Dks){
  ##############################################################################
  # Dimensions
  # Number of Observations
  N  <- dim(y)[1]
  KK <- dim(y)[2]
  ##############################################################################
  
  ##############################################################################
  # Auxiliary Variables
  # Transforms Back to regular variables
  t2 <- exp( ga) * s2
  l2 <- exp(-ph)
  # Covariances
  Ckk  <- t2 * exp(- Dkk / l2)
  Cks  <- t2 * exp(- Dks / l2)
  K    <- Ckk + Cks %*% t(Cks) * KK / s2
  eK   <- eigen(K)
  Cy   <- Cks %*% yb
  eCkk <- eigen(Ckk)
  # Loglikelihood
  l <-     - sum(log(eK$values))/2 + sum(log(eCkk$values))/2
  return(l)
  l <- l + t(Cy) %*% eK$vectors %*% diag(1 / eK$values) %*% t(eK$vectors) %*% Cy / (2 * s2^2)
  l <- l + ga - ph
  ##############################################################################
  
  ##############################################################################
  # Returns
  return(l)
  ##############################################################################
}