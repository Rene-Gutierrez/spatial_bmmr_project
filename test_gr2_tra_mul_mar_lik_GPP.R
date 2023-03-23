################################################################################
# Test Gradient Marginal Likelihood Gaussian Predictive Process (Tra and Mul)
################################################################################

# Libraries
source("./gra_tr2_mul_mar_lik_GPP.R")

# Settings
h   <- 0.00000001
# Real Scale
ga  <-   log(t2 / s2)
ph  <-   1  / l2
gah <- ga + h
phh <- ph + h
yb  <- rowMeans(y)

# Regular Expressions
Ckk  <- exp(- ph * Dkk)
Cks  <- exp(- ph * Dks)
A    <- s2 * (exp(ga) * t(Cks) %*% solve(Ckk, Cks) + diag(nrow = N) / K)
eig  <- eigen(A)

# # Test ga
# Ckkh <- exp(- ph * Dkk)
# Cksh <- exp(- ph * Dks)
# Ah   <- s2 * (exp(gah) * t(Cksh) %*% solve(Ckkh, Cksh) + diag(nrow = N) / K)
# eigh <- eigen(Ah)
# out2 <- ((- sum(log(eigh$values)) / 2) - (- sum(log(eig$values)) / 2)) / h
# out3 <- ((- t(yb) %*% solve(Ah, yb) / 2) - (- t(yb) %*% solve(A, yb) / 2)) / h
# out4 <- out2 + out3 + 1

# Test ph
Ckkh <- exp(- phh * Dkk)
Cksh <- exp(- phh * Dks)
Ah   <- s2 * (exp(ga) * t(Cksh) %*% solve(Ckkh, Cksh) + diag(nrow = N) / K)
eigh <- eigen(Ah)
out5 <- gra_tr2_mul_mar_lik_GPP(s2  = s2,
                                ga  = ga,
                                ph  = ph,
                                yb  = yb,
                                Dkk = Dkk,
                                Dks = Dks,
                                KK  = K)
out6 <- ((- sum(log(eigh$values)) / 2) - (- sum(log(eig$values)) / 2)) / h
out7 <- ((- t(yb) %*% solve(Ah, yb) / 2) - (- t(yb) %*% solve(A, yb) / 2)) / h
out8 <- out6 + out7


tim <- Sys.time()
sg  <- matrix(data = NA, nrow = I, ncol = 2)
sph <- exp(seq(log(0.1), log(10), length.out = 100))
for(i in 1:100){
  g <- gra_tr2_mul_mar_lik_GPP(s2  = s2,
                               ga  = ga,
                               ph  = sph[i],
                               yb  = yb,
                               Dkk = Dkk,
                               Dks = Dks,
                               KK  = K)
  sg[i, ] <- g
}
tim <- Sys.time() - tim
print(tim)
