################################################################################
# MALA Step (Transformed Covariance Parameters)
################################################################################

mala_ste_tra_gp <- function(ls2,
                            lt2,
                            ll2,
                            y,
                            DM,
                            h = 1){
  ##############################################################################
  # MALA Step
  # Current Observation
  xc  <- c(ls2, lt2, ll2)
  # Current Gradient
  gxc <- tra_gra_mar_lik_GP(ls2 = ls2,
                            lt2 = lt2,
                            ll2 = ll2,
                            y  = y,
                            DM = DM)
  # Proposal
  xp  <- xc + h * gxc + rnorm(n = 3, sd = sqrt(2 * h))
  # Proposal Gradient
  gxp <- tra_gra_mar_lik_GP(ls2 = xp[1],
                            lt2 = xp[2],
                            ll2 = xp[3],
                            y  = y,
                            DM = DM)
  # Loglikelihood Proposal
  lxp <- tra_mar_lik_GP(ls2 = xp[1],
                        lt2 = xp[2],
                        ll2 = xp[3],
                        y  = y,
                        DM = DM)
  # Loglikelihood Current
  lxc <- tra_mar_lik_GP(ls2 = ls2,
                        lt2 = lt2,
                        ll2 = ll2,
                        y  = y,
                        DM = DM)
  # Acceptance Rate
  aR <- (lxp - lxc) - (sum((xc - xp - h * gxp)^2) - sum((xp - xc - h * gxc)^2)) / (4 * h)
  aR <- exp(aR)
  # Checks Acceptance
  if(aR > runif(n = 1)){
    xn <- xp
  } else {
    xn <- xc
  }
  # Return
  return(xn)
}