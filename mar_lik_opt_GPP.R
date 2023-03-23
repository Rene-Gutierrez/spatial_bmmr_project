################################################################################
# Mrginal Likelihood Optimization GPP
################################################################################

h  <- 0.001
S  <- 300
x  <- c(2, 0.1, 2)
sx <- matrix(data = NA, nrow = S, ncol = 3)

tim <- Sys.time()
for(i in 1:S){
  # print(i)
  g <- gra_mar_lik_GPP(s2  = x[1],
                       t2  = x[2],
                       l2  = x[3],
                       y   = y,
                       Dkk = Dkk,
                       Dks = Dks,
                       I   = 10)
  x <- x + h * g
  sx[i, ] <- x
}
tim <- Sys.time() - tim

plot(sx[, 1], type = 'l', ylim = c(min(s2, sx[, 1]), max(s2, sx[, 1])))
abline(h = s2, col = 'red')
plot(sx[, 2], type = 'l', ylim = c(min(t2, sx[, 2]), max(t2, sx[, 2])))
abline(h = t2, col = 'red')
plot(sx[, 3], type = 'l', ylim = c(min(l2, sx[, 3]), max(l2, sx[, 3])))
abline(h = l2, col = 'red')

print(tim)
