################################################################################
# Predictive Gaussian Process Sampler for w tilde
################################################################################

w_tilde_sampler <- function(s2,
                            ga,
                            ph,
                            Dkk,
                            Dks,
                            yb,
                            K){
  # Auxiliary Variables
  M    <- dim(Dkk)[1]
  t2   <- ga * s2
  Cks  <- exp(- ph * Dks)
  Ckk  <- exp(- ph * Dkk)
  eCkk <- eigen(Ckk)
  iC   <- eCkk$vectors %*% (t(eCkk$vectors) / eCkk$values)
  CC   <- iC %*% Cks
  Sig  <- (CC %*% t(CC) * (K / s2) + iC / t2)
  eSig <- eigen(Sig)
  iSig <- eSig$vectors %*% (t(eSig$vectors) / eSig$values)
  # print(dim(iSig))
  # print(dim(iC))
  # print(dim(Cks))
  # print(length(yb))
  wb   <- iSig %*% iC %*% Cks %*% yb * (K / s2)
  # Samples w star
  w  <- rnorm(n = M, mean = 0, sd = 1)
  w  <- eSig$vectors %*% (t(eSig$vectors) * sqrt(1 / eSig$values)) %*% w
  w  <- w + wb
  wt <- t(CC) %*% w
  return(wt)
}