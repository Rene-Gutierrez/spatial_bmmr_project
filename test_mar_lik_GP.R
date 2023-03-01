################################################################################
# Marginal Loglikelihood Test
################################################################################

out1 <- mar_lik_GP(s2 = s2,
                   t2 = t2,
                   l2 = l2,
                   y  = y,
                   DM = DM)

out2 <- mvtnorm::dmvnorm(x     = y,
                         mean  = rep(0, N),
                         sigma = t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N),
                         log   = TRUE)
