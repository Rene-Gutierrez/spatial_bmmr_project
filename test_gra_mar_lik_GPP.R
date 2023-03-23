################################################################################
# Test Gradient Marginal Likelihood Gaussian Predictive Process
################################################################################

# Settings
h <- 0.00000001

# Regular Expressions
Ckk  <- t2 * exp(- Dkk / l2)
Cks  <- t2 * exp(- Dks / l2)
C    <- t(Cks) %*% solve(Ckk, Cks) + s2 * diag(nrow = N)
eig  <- eigen(C)

# Tests s2
Ckkh <- t2 * exp(- Dkk / l2)
Cksh <- t2 * exp(- Dks / l2)
Ch   <- t(Cksh) %*% solve(Ckkh, Cksh) + (s2 + h) * diag(nrow = N)
eigh <- eigen(Ch)
out1 <- gra_mar_lik_GPP(s2  = s2,
                        t2  = t2,
                        l2  = l2,
                        y   = y,
                        Dkk = Dkk,
                        Dks = Dks)
out2 <- ((- sum(log(eigh$values)) / 2) - (- sum(log(eig$values)) / 2)) / h
out3 <- ((- t(y) %*% solve(Ch, y) / 2) - (- t(y) %*% solve(C, y) / 2)) / h
out4 <- mar_lik_GPP(s2  = s2 + h,
                    t2  = t2,
                    l2  = l2,
                    y   = y,
                    Dkk = Dkk,
                    Dks = Dks) - mar_lik_GPP(s2  = s2,
                                             t2  = t2,
                                             l2  = l2,
                                             y   = y,
                                             Dkk = Dkk,
                                             Dks = Dks)
out4 <- out4 / h

# Tests t2
Ckkh <- (t2 + h) * exp(- Dkk / l2)
Cksh <- (t2 + h) * exp(- Dks / l2)
Ch   <- t(Cksh) %*% solve(Ckkh, Cksh) + s2 * diag(nrow = N)
eigh <- eigen(Ch)
out1 <- gra_mar_lik_GPP(s2  = s2,
                        t2  = t2,
                        l2  = l2,
                        y   = y,
                        Dkk = Dkk,
                        Dks = Dks)
out2 <- ((- sum(log(eigh$values)) / 2) - (- sum(log(eig$values)) / 2)) / h
out3 <- ((- t(y) %*% solve(Ch, y) / 2) - (- t(y) %*% solve(C, y) / 2)) / h
out7 <- mar_lik_GPP(s2  = s2,
                    t2  = t2 + h,
                    l2  = l2,
                    y   = y,
                    Dkk = Dkk,
                    Dks = Dks) - mar_lik_GPP(s2  = s2,
                                             t2  = t2,
                                             l2  = l2,
                                             y   = y,
                                             Dkk = Dkk,
                                             Dks = Dks)
out7 <- out7 / h

# Tests l
Ckkh <- t2 * exp(- Dkk / (l2 + h))
Cksh <- t2 * exp(- Dks / (l2 + h))
Ch   <- t(Cksh) %*% solve(Ckkh, Cksh) + s2 * diag(nrow = N)
eigh <- eigen(Ch)
out1 <- gra_mar_lik_GPP(s2  = s2,
                        t2  = t2,
                        l2  = l2,
                        y   = y,
                        Dkk = Dkk,
                        Dks = Dks)
out2 <- ((- sum(log(eigh$values)) / 2) - (- sum(log(eig$values)) / 2)) / h
out3 <- ((- t(y) %*% solve(Ch, y) / 2) - (- t(y) %*% solve(C, y) / 2)) / h
out10 <- mar_lik_GPP(s2  = s2,
                    t2  = t2,
                    l2  = l2 + h,
                    y   = y,
                    Dkk = Dkk,
                    Dks = Dks) - mar_lik_GPP(s2  = s2,
                                             t2  = t2,
                                             l2  = l2,
                                             y   = y,
                                             Dkk = Dkk,
                                             Dks = Dks)
out10 <- out10 / h
