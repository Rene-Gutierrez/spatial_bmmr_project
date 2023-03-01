################################################################################
# Test Gradient Marginal Likelihood Gaussian Process
################################################################################
h    <- 0.00000001
out1 <- gra_mar_lik_GP(s2 = s2,
                       t2 = t2,
                       l2 = l2,
                       y  = y,
                       DM = DM)
out2 <- mar_lik_GP(s2 = s2,
                   t2 = t2,
                   l2 = l2,
                   y  = y,
                   DM = DM) - 
  mar_lik_GP(s2 = s2 + h,
             t2 = t2,
             l2 = l2,
             y  = y,
             DM = DM)
out2 <- out2 / h

out3 <- - log(det((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)))) / 2 + 
  log(det((t2 * exp(- DM / (2 * l2)) + (s2 + h) * diag(nrow = N)))) / 2
out3 <- out3 / h

out4 <- - t(y) %*% solve((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)), y) / 2 -
  (- t(y) %*% solve((t2 * exp(- DM / (2 * l2)) + (s2 + h) * diag(nrow = N)), y) / 2)
out4 <- out4 / h

C <- t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)
svC <- svd(C)

out5 <- mvtnorm::dmvnorm(x     = y,
                         mean  = rep(0, N),
                         sigma = t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N),
                         log   = TRUE) - mvtnorm::dmvnorm(x     = y,
                                                          mean  = rep(0, N),
                                                          sigma = t2 * exp(- DM / (2 * l2)) + (s2 + h) * diag(nrow = N),
                                                          log   = TRUE)
out5 <- out5 / h

out6 <- - log(det(solve((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N))))) / 2 - 
  t(y) %*% solve((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)), y) / 2

out7 <- - log(det((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)))) / 2 -
  t(y) %*% solve((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)), y) / 2 - N * log(2 * pi) / 2


out8 <- - log(det((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)))) / 2 + 
  log(det(((t2 + h) * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)))) / 2
out8 <- out8 / h


out9 <- - t(y) %*% solve((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)), y) / 2 -
  (- t(y) %*% solve(((t2 + h) * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)), y) / 2)
out9 <- out9 / h

out10 <- mar_lik_GP(s2 = s2,
                   t2 = t2,
                   l2 = l2,
                   y  = y,
                   DM = DM) - 
  mar_lik_GP(s2 = s2,
             t2 = t2 + h,
             l2 = l2,
             y  = y,
             DM = DM)
out10 <- out10 / h

out11 <- - log(det((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)))) / 2 + 
  log(det((t2 * exp(- DM / (2 * (l2 + h))) + s2 * diag(nrow = N)))) / 2
out11 <- out11 / h

out12 <- - t(y) %*% solve((t2 * exp(- DM / (2 * l2)) + s2 * diag(nrow = N)), y) / 2 -
  (- t(y) %*% solve((t2 * exp(- DM / (2 * (l2 + h))) + s2 * diag(nrow = N)), y) / 2)
out12 <- out12 / h

out13 <- mar_lik_GP(s2 = s2,
                    t2 = t2,
                    l2 = l2,
                    y  = y,
                    DM = DM) - 
  mar_lik_GP(s2 = s2,
             t2 = t2,
             l2 = l2 + h,
             y  = y,
             DM = DM)
out13 <- out13 / h

