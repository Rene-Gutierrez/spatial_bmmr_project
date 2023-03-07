################################################################################
# Gradient Marginal Likelihood Gaussian Predictive Process
################################################################################

gra_mar_lik_GPP <- function(s2,
                            t2,
                            l2,
                            y,
                            Dkk,
                            Dks){
  ##############################################################################
  # Dimensions
  # Number of Observations
  N   <- length(y)
  ##############################################################################
  
  ##############################################################################
  # Auxiliary Variables
  # Covariances
  Ckk  <- t2 * exp(- Dkk / l2)
  Cks  <- t2 * exp(- Dks / l2)
  K    <- Ckk + Cks %*% t(Cks) / s2
  eK   <- eigen(K)
  Cy   <- Cks %*% y
  iK   <- eK$vectors %*% diag(1 / eK$values) %*% t(eK$vectors)
  iC   <- diag(nrow = N) / s2 - t(Cks) %*% iK %*% Cks / (s2^2)
  Cy   <- iC %*% y
  CCs  <- solve(Ckk, Cks)
  CCC  <- t(Cks) %*% solve(Ckk, Cks)
  iCkk <- solve(Ckk)
  # Auxiliary Derivatives
  dCks  <-   Cks * Dks / (l2^2)
  dCkk  <-   Ckk * (Dkk / (l2^2))
  diCkk <- - iCkk %*% dCkk %*% iCkk
  dC    <-   t(Cks) %*% (iCkk %*% dCks + diCkk %*% Cks) + t(dCks) %*% iCkk %*% Cks
  # Derivative s2
  gs2 <- - N / s2 / 2
  gs2 <- gs2 + sum(diag(t(Cks) %*% iK %*% Cks)) / (2 * s2^2)
  gs2 <- gs2 + t(Cy) %*% Cy / 2
  # Derivative t2
  gt2 <- - sum(diag(iC %*% t(Cks) %*% solve(Ckk, Cks))) / t2 / 2
  gt2 <- gt2 + t(Cy) %*% CCC %*% Cy / t2 / 2
  # Derivative l2
  gl2 <-     - sum(diag(iC %*% dC)) / 2
  gl2 <- gl2 + y %*% iC %*% dC %*% iC %*% y / 2
  # Returns
  return(c(gs2, gt2, gl2))
}