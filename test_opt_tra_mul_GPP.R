################################################################################
# Test Optimization Step (Tra and Mul)
################################################################################
I   <- 300
h   <- 0.01
ga  <- t2 / s2
ph  <- 1  / l2
x   <- c(1, 1)
sx  <- matrix(data = NA, nrow = I, ncol = 2)
sg  <- matrix(data = NA, nrow = I, ncol = 2)
tim <- Sys.time()
for(i in 1:I){
  g <- gra_tra_mul_mar_lik_GPP(s2  = s2,
                               ga  = x[1],
                               ph  = x[2],
                               yb  = yb,
                               Dkk = Dkk,
                               Dks = Dks,
                               KK  = K)
  x <- x + h * g
  sg[i, ] <- g
  sx[i, ] <- x
}
tim <- Sys.time() - tim
print(tim)

plot(sx[,1], type = 'l')
abline(h = ga, col = 'red')
plot(sx[,2], type = 'l')
abline(h = ph, col = 'red')
