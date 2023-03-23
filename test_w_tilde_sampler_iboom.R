source("./w_tilde_sampler.R")

p  <- 7

yb <- colMeans(G[,C[,p],p])
s2 <- 1
ga <- 0.5
ph <- 1

out <- w_tilde_sampler(s2  = s2,
                       ga  = ga,
                       ph  = 1,
                       Dkk = Dkk[[p]],
                       Dks = Dks[[p]],
                       yb  = yb,
                       K   = dim(G[,C[,p],p])[1])

ymax <- max(yb, out$wt)
ymin <- min(yb, out$wt)

plot(yb, type = 'l', ylim = c(ymin, ymax))
par(new = TRUE)
plot(out$wt, type = 'l', ylim = c(ymin, ymax), col = 'red')
