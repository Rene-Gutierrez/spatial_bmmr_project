# Auxiliary Variables
out <- w_tilde_sampler(s2  = s2,
                       ga  = t2 * s2,
                       ph  = 1 / l2,
                       Dkk = Dkk,
                       Dks = Dks)

# Plotting Dimensions
ymax <- max(   3)
ymin <- min(  -3)
# Plotting Order
o    <- order(S[,1])
# Plots
plot(x    = S[,1][o],
     y    = w[o],
     ylim = c(ymin, ymax),
     type = 'l')
for(k in 1:K){
  par(new=TRUE)
  plot(x    = S[,1][o],
       y    = y[o, k],
       ylim = c(ymin, ymax),
       col  = rgb(0, 0, 1, 0.1),
       type = 'l')
}
par(new=TRUE)
plot(x    = S[,1][o],
     y    = w[o],
     ylim = c(ymin, ymax),
     type = 'l')
par(new=TRUE)
plot(x    = S[,1][o],
     y    = yb[o],
     ylim = c(ymin, ymax),
     type = 'l')
par(new=TRUE)
plot(x    = Sk,
     y    = out$w,
     ylim = c(ymin, ymax),
     type = 'l',
     col  = 'green',
     lwd  = 2)
par(new=TRUE)
plot(x    = Sk,
     y    = out$wb,
     ylim = c(ymin, ymax),
     type = 'l',
     col  = 'red',
     lwd  = 2)
# par(new=TRUE)
# plot(x    = Sk,
#      y    = out$wb + out$w,
#      ylim = c(ymin, ymax),
#      type = 'l',
#      col  = 'purple',
#      lwd  = 2)
# plot(x    = S[,1][o],
#      y    = colMeans(swt)[o],
#      ylim = c(ymin, ymax),
#      type = 'l',
#      lwd  = 2,
#      col  = 'red')
# par(new=TRUE)
# plot(x    = S[,1][o],
#      y    = swt[1,][o],
#      ylim = c(ymin, ymax),
#      type = 'l',
#      lwd  = 2,
#      col  = 'green')
