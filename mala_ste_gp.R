################################################################################
# MALA Step
################################################################################

mala_ste_gp <- function(s2,
                        t2,
                        l2,
                        y,
                        DM,
                        h = 1){
  ##############################################################################
  # MALA Step
  # Current Observation
  xc  <- c(s2, t2, l2)
  # Current Gradient
  gxc <- gra_mar_lik_GP(s2 = s2,
                        t2 = t2,
                        l2 = l2,
                        y  = y,
                        DM = DM)
  # Proposal
  xp  <- xc + h * gxc + rnorm(n = 3, sd = sqrt(2 * h))
  # Proposal Gradient
  gxp <- gra_mar_lik_GP(s2 = xp[1],
                        t2 = xp[2],
                        l2 = xp[3],
                        y  = y,
                        DM = DM)
  # Loglikelihood Proposal
  lxp <- mar_lik_GP(s2 = xp[1],
                    t2 = xp[2],
                    l2 = xp[3],
                    y  = y,
                    DM = DM)
  # Loglikelihood Current
  lxc <- mar_lik_GP(s2 = s2,
                    t2 = t2,
                    l2 = l2,
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