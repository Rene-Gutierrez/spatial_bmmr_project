################################################################################
# Gradient Marginal Likelihood Gaussian Predictive Process (Mul and Tra)
################################################################################

gra_tr2_mul_mar_lik_GPP <- function(s2,
                                    ga,
                                    ph,
                                    yb,
                                    Dkk,
                                    Dks,
                                    KK){
  ##############################################################################
  # Dimensions
  # Number of Observations
  N <- length(yb)
  ##############################################################################
  
  ##############################################################################
  # Auxiliary Variables
  # Covariances
  Ckk  <- exp(- ph * Dkk)
  Cks  <- exp(- ph * Dks)
  K    <- Ckk / exp(ga) + Cks %*% t(Cks) * KK
  iC   <- diag(nrow = N) * KK - t(Cks) %*% solve(K, Cks) * KK^2
  Cy   <- iC %*% yb
  CCC  <- t(Cks) %*% solve(Ckk, Cks)
  iCkk <- solve(Ckk)
  # Auxiliary Derivatives
  dCks  <- - Cks * (Dks) 
  dCkk  <- - Ckk * (Dkk)
  diCkk <- - iCkk %*% dCkk %*% iCkk
  dC    <-   t(Cks) %*% (iCkk %*% dCks + diCkk %*% Cks) + t(dCks) %*% iCkk %*% Cks
  dC    <- exp(ga) * dC
  # Derivative ga
  gga <-     - sum(iC * CCC) * exp(ga) / 2
  gga <- gga + t(Cy) %*% CCC %*% Cy * exp(ga) / 2 / s2
  gga <- gga + 1
  # Derivative l2
  gph <- - sum(iC * dC) / 2
  gph <- gph + t(Cy) %*% dC %*% Cy / 2 / s2
  gph <- gph
  return(c(gga, gph))
}