################################################################################
### Saptial Data Run
################################################################################

tim_set <- Sys.time()
################################################################################
# Libraries and Functions
source("./group_iterator.R")
source("./spa_iboom_iterator.R")
source("./spa_iboom_sampler.R")
source("./w_tilde_sampler.R")
library(mvtnorm)
################################################################################

################################################################################
# Set-Up
# Seed
set.seed(25052022)
# SPN Region Selection
spn <- c(0, 1, 2, 3)
tag <- c()
for(i in spn){
  tag <- paste0(tag, i)
}
# Outcome Selection
var <- 48
# Number of Knots
K   <- 125
################################################################################

################################################################################
# Data Import
A      <- readRDS(file = paste0("./dat/rea_A_", tag, ".rds"))
G      <- readRDS(file = paste0("./dat/rea_G_", tag, ".rds"))
y      <- readRDS(file = "./dat/rea_y.rds")
C      <- readRDS(file = paste0("./dat/rea_C_", tag, ".rds"))
L      <- readRDS(file = paste0("./dat/rea_L_", tag, ".rds"))
Dkk    <- readRDS(file = paste0("./dat/Dkk_", K, "_", tag, ".rds"))
Dks    <- readRDS(file = paste0("./dat/Dks_", K, "_", tag, ".rds"))
reg_id <- readRDS(file = paste0("./dat/regID_", tag, ".rds"))
################################################################################

################################################################################
# Data Wrangling
# Covariate Selection
X   <- y[, c(5, 12, 13)]
# Speech Outcome
y   <- y[, var]
# Removes NA
selX <- rowSums(is.nan(as.matrix(X))) == 0
sel  <- !is.na(y)
sel  <- sel & selX
y    <- y[sel]
X    <- X[sel,]
A    <- A[sel, ,]
G    <- G[sel, ,]
# Objects Dimensions
N  <- length(y)
P  <- dim(G)[3]
mV <- dim(G)[2]

# Normalizes
# Centers
mA <- apply(X = A, MARGIN = c(2, 3), FUN = mean)
# mG <- apply(X = G, MARGIN = c(2, 3), FUN = mean)
A  <- A - array(data = mA %x% rep(1, N), dim = c(N, P,  P))
# G  <- G - array(data = mG %x% rep(1, N), dim = c(N, mV, P))
G  <- G - mean(G, na.rm = TRUE)
# Normalizes
y  <- (y - mean(y)) / sd(y)
# Vectorizes A and G
Zv <- matrix(data = NA, nrow = N, ncol = sum(C) + P * (P - 1) / 2)
for(i in 1:N){
  Zv[i, ] <- c(A[i,,][upper.tri(A[i,,])], G[i,,][C])
}
# Generates Covariates
X[, 1] <- (X[, 1] == 1) + 0
X[, 2] <- X[, 2] - mean(X[, 2])
X[, 3] <- X[, 3] - mean(X[, 3])
X      <- as.matrix(X)
XG     <- kronecker(X = rep(1, mV), Y = X)
saveRDS(object = X, file = paste0("./dat/rea_X_", var, ".rds"))
saveRDS(object = y, file = paste0("./dat/rea_y_", var, ".rds"))
saveRDS(object = A, file = paste0("./dat/rea_A_", var, "_", tag, ".rds"))
saveRDS(object = G, file = paste0("./dat/rea_G_", var, "_", tag, ".rds"))
################################################################################
tim_set <- Sys.time() - tim_set
print(tim_set)

tim_run <- Sys.time()
################################################################################
# Run Set-Up
# Number of Samples (After Burn-In)
nmcmc  <- 300
# Number of Burn-In Samples
nburn  <- 100
# Number of Batches to Divide the Samples
nbatch <- 10
# Size of the Batch
batSiz <- nmcmc / nbatch
################################################################################

################################################################################
# Run
# First Run
print(paste0("Number of Batches = ", nbatch))
# First Run
print(paste0("Batch Number = 1"))
out_rea <- spa_iboom_sampler(y      = y,
                             X      = X,
                             XG     = XG,
                             Zv     = Zv,
                             G      = G,
                             A      = A,
                             C      = C,
                             Dkk    = Dkk,
                             Dks    = Dks,
                             nmcmc  = batSiz,
                             burnin = nburn,
                             onlyg  = TRUE,
                             state  = list())
print("")
# Saves Samples
saveRDS(object = nbatch,         file = paste0("./out/nbatch", var, "_", tag, ".rds"))
saveRDS(object = batSiz,         file = paste0("./out/batSiz", var, "_", tag, ".rds"))
saveRDS(object = out_rea$sam$g,  file = paste0("./out/g_",     var, "_", tag, "_", "batch=", "1", ".rds"))
saveRDS(object = out_rea$sam$B,  file = paste0("./out/B_",     var, "_", tag, "_", "batch=", "1", ".rds"))
saveRDS(object = out_rea$sam$DA, file = paste0("./out/DA_",    var, "_", tag, "_", "batch=", "1", ".rds"))
saveRDS(object = out_rea$sam$DG, file = paste0("./out/DG_",    var, "_", tag, "_", "batch=", "1", ".rds"))
saveRDS(object = out_rea$sam$wG, file = paste0("./out/wG_",    var, "_", tag, "_", "batch=", "1", ".rds"))
saveRDS(object = out_rea$sam$s2, file = paste0("./out/s2_",    var, "_", tag, "_", "batch=", "1", ".rds"))
# Rest of Runs
for(batch in 1:nbatch){
  print(paste0("Batch Number = ", batch))
  # First Run
  out_rea <- spa_iboom_sampler(y      = y,
                               X      = X,
                               XG     = XG,
                               Zv     = Zv,
                               G      = G,
                               A      = A,
                               C      = C,
                               Dkk    = Dkk,
                               Dks    = Dks,
                               nmcmc  = batSiz,
                               burnin = 0,
                               onlyg  = TRUE,
                               state  = out_rea$state)
  print("")

  saveRDS(object = out_rea$sam$g,  file = paste0("./out/g_",     var, "_", tag, "_", "batch=", batch, ".rds"))
  saveRDS(object = out_rea$sam$B,  file = paste0("./out/B_",     var, "_", tag, "_", "batch=", batch, ".rds"))
  saveRDS(object = out_rea$sam$DA, file = paste0("./out/DA_",    var, "_", tag, "_", "batch=", batch, ".rds"))
  saveRDS(object = out_rea$sam$DG, file = paste0("./out/DG_",    var, "_", tag, "_", "batch=", batch, ".rds"))
  saveRDS(object = out_rea$sam$wG, file = paste0("./out/wG_",    var, "_", tag, "_", "batch=", batch, ".rds"))
  saveRDS(object = out_rea$sam$s2, file = paste0("./out/s2_",    var, "_", tag, "_", "batch=", batch, ".rds"))
}
################################################################################
tim_run <- Sys.time() - tim_run
print(tim_run)