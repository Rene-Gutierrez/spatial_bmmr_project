################################################################################
# Test Transformed MALA Step
################################################################################
S   <- 500

out <- mala_ste_tra_gp(ls2 = 0,
                       lt2 = 0,
                       ll2 = 0,
                       y   = y,
                       DM  = DM,
                       h   = 1)

sP <- matrix(data = NA, nrow = S, ncol = 3)
for(i in 1:S){
  out <- mala_ste_tra_gp(ls2 = out[1],
                         lt2 = out[2],
                         ll2 = out[3],
                         y  = y,
                         DM = DM,
                         h  = 0.1)
  sP[i, ] <- out
}

# Plots
plot(sP[, 1],
     type = 'l',
     ylim = c(min(log(s2), sP[, 1]), max(log(s2), sP[, 1])))
abline(h = log(s2), col = 'red')
plot(exp(sP[, 1]),
     type = 'l',
     ylim = c(min(s2, exp(sP[, 1])), max(s2, exp(sP[, 1]))))
abline(h = s2, col = 'red')
plot(exp(sP[, 2]),
     type = 'l',
     ylim = c(min(t2, exp(sP[, 2])), max(t2, exp(sP[, 2]))))
abline(h = t2, col = 'red')
plot(exp(sP[, 3]),
     type = 'l',
     ylim = c(min(l2, exp(sP[, 3])), max(l2, exp(sP[, 3]))))
abline(h = l2, col = 'red')
