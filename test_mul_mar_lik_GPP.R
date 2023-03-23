################################################################################
# Test Gaussian Predictive Process Marginal Likelihood
################################################################################

source("./mul_mar_lik_GPP.R")
source("./mar_lik_GPP.R")

Ckk  <- t2 * exp(- Dkk / l2)
Cks  <- t2 * exp(- Dks / l2)

eig  <- eigen(t(Cks) %*% solve(Ckk, Cks) + s2 * diag(nrow = N))
out1 <- - sum(log(eig$values)) / 2
out2 <- - t(y) %*% eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors) %*% y / 2
tim1 <- Sys.time()
out3 <- - t(y) %*% solve(t(Cks) %*% solve(Ckk, Cks) + s2 * diag(nrow = N), y) / 2
tim1 <- Sys.time() - tim1

tim2 <- Sys.time()
out4 <- mul_mar_lik_GPP(s2  = s2,
                        t2  = t2,
                        l2  = l2,
                        y   = y,
                        Dkk = Dkk,
                        Dks = Dks)
tim2 <- Sys.time() - tim2

tim3 <- Sys.time()
out5 <- 0
for(k in 1:K){
  out5 <- out5 + mar_lik_GPP(s2  = s2,
                             t2  = t2,
                             l2  = l2,
                             y   = y[,k],
                             Dkk = Dkk,
                             Dks = Dks)
}
tim3 <- Sys.time() - tim3