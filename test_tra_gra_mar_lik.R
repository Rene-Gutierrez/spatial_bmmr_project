################################################################################
# Test Transformed Gradient Marginal Likelihood
################################################################################

# Library
source("./gaussian_process_simulation.R")

ls2  <- log(s2)
lt2  <- log(t2)
ll2  <- log(l2)
h    <- 0.00000001

# Auxilary Variables
K    <- exp(lt2) * exp(- DM / (2 * exp(ll2)))
C    <- K + exp(ls2) * diag(nrow = N)

# Test for ls2
Kh   <- exp(lt2) * exp(- DM / (2 * exp(ll2)))
Ch   <- K + exp(ls2 + h) * diag(nrow = N)

# Test Derivative of C (ls2)
dC_app <- (Ch - C) / h
dC_exa <- exp(ls2) * diag(nrow = N)

# Test Derivative of log|C| (ls2)
dlogdetC_app <- (log(det(Ch)) - log(det(C))) / h
dlogdetC_exa <- exp(ls2) * sum(diag(solve(C)))

# Test Derivative of ls2 
out1 <- tra_gra_mar_lik_GP(ls2 = ls2,
                           lt2 = lt2,
                           ll2 = ll2,
                           y   = y,
                           DM  = DM)
out2 <- (- log(det(Ch)) / 2 - (- log(det(C))) / 2) / h
out3 <- ((- t(y) %*% solve(Ch) %*% y / 2) - (- t(y) %*% solve(C) %*% y / 2)) / h
out4 <- (tra_mar_lik_GP(ls2 = ls2 + h,
                        lt2 = lt2,
                        ll2 = ll2,
                        y   = y,
                        DM  = DM) - tra_mar_lik_GP(ls2 = ls2,
                                                   lt2 = lt2,
                                                   ll2 = ll2,
                                                   y   = y,
                                                   DM  = DM)) / h

# Test for lt2
Kh   <- exp(lt2 + h) * exp(- DM / (2 * exp(ll2)))
Ch   <- Kh + exp(ls2) * diag(nrow = N)

out5 <- (- log(det(Ch)) / 2 - (- log(det(C))) / 2) / h
out6 <- ((- t(y) %*% solve(Ch) %*% y / 2) - (- t(y) %*% solve(C) %*% y / 2)) / h
out7 <- (tra_mar_lik_GP(ls2 = ls2,
                        lt2 = lt2 + h,
                        ll2 = ll2,
                        y   = y,
                        DM  = DM) - tra_mar_lik_GP(ls2 = ls2,
                                                   lt2 = lt2,
                                                   ll2 = ll2,
                                                   y   = y,
                                                   DM  = DM)) / h

# Test for ll2
Kh   <- exp(lt2) * exp(- DM / (2 * exp(ll2 + h)))
Ch   <- Kh + exp(ls2) * diag(nrow = N)

out8  <- (- log(det(Ch)) / 2 - (- log(det(C))) / 2) / h
out9  <- ((- t(y) %*% solve(Ch) %*% y / 2) - (- t(y) %*% solve(C) %*% y / 2)) / h
out10 <- (tra_mar_lik_GP(ls2 = ls2,
                         lt2 = lt2,
                         ll2 = ll2 + h,
                         y   = y,
                         DM  = DM) - tra_mar_lik_GP(ls2 = ls2,
                                                    lt2 = lt2,
                                                    ll2 = ll2,
                                                    y   = y,
                                                    DM  = DM)) / h