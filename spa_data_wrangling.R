### SPN Selection
tag <- c()
for(i in spn){
  tag <- paste0(tag, i)
}
reg <- readRDS(file = paste0("regID_", tag, ".rds"))

### Outcome
#var <- 42

### Data Import
nbatch <- readRDS(file = paste0("./out/nbatch", var, "_", tag, ".rds"))
batSiz <- readRDS(file = paste0("./out/batSiz", var, "_", tag, ".rds"))
C      <- readRDS(file = paste0("./dat/rea_C_", tag, ".rds"))
# Problem Dimensions
P  <- dim(C)[2]
mV <- dim(C)[1]

# Imports B
# First Batch
B    <- readRDS(file = paste0("./out/B_",     var, "_", tag, "_", "batch=", "1", ".rds"))
# Rest of Batches
for(batch in 2:nbatch){
  B    <- rbind(B,
                readRDS(file = paste0("./out/B_",     var, "_", tag, "_", "batch=", batch, ".rds")))
}
Bq   <- matrix(data = NA, nrow = dim(B)[2], ncol = 3)
Bq[, 1] <- apply(X = B, MARGIN = 2, FUN = quantile, probs = 0.025)
Bq[, 2] <- apply(X = B, MARGIN = 2, FUN = quantile, probs = 0.500)
Bq[, 3] <- apply(X = B, MARGIN = 2, FUN = quantile, probs = 0.975)
remove(B)
gc()

# Theta
Theta <- array(data = 0, dim = c(3, P, P))
for(i in 1:3){
  Theta[i,,][upper.tri(Theta[i,,])] <- Bq[1:(P * (P - 1) / 2), i]
  Theta[i,,]                        <- Theta[i,,] + t(Theta[i,,])
  diag(Theta[i,,]) <- NA
}
# Beta
Beta <- array(data = NA, dim = c(3, mV, P))
for(i in 1:3){
  Beta[i,,][C] <- Bq[-(1:(P * (P - 1) / 2)), i]
}
remove(Bq)
gc()
# Saves the Values
saveRDS(object = Theta,  file = paste0("./out/Theta_q_", var, "_", tag, ".rds"))
saveRDS(object = Beta,   file = paste0("./out/Beta_q_",  var, "_", tag, ".rds"))

# Imports DA and DG
M   <- dim(readRDS(file = paste0("./out/DA_",     var, "_", tag, "_", "batch=", "1", ".rds")))[3]
DAq <- array(data = 0, dim = c(3, P, M))
DGq <- array(data = 0, dim = c(3, P, M))
for(m in 1:M){
  print(m)
  # First Batch
  DA <- readRDS(file = paste0("./out/DA_",     var, "_", tag, "_", "batch=", "1", ".rds"))[,,m]
  # Rest of Batches
  for(batch in 2:nbatch){
    DA   <- rbind(DA,
                  readRDS(file = paste0("./out/DA_",     var, "_", tag, "_", "batch=", batch, ".rds"))[,,m])
  }
  DAq[1, , m] <- apply(X = DA, MARGIN = 2, FUN = quantile, probs = 0.025)
  DAq[2, , m] <- apply(X = DA, MARGIN = 2, FUN = quantile, probs = 0.500)
  DAq[3, , m] <- apply(X = DA, MARGIN = 2, FUN = quantile, probs = 0.975)
  remove(DA)
  gc()
  
  # First Batch
  DG <- readRDS(file = paste0("./out/DG_",     var, "_", tag, "_", "batch=", "1", ".rds"))[,,m]
  # Rest of Batches
  for(batch in 2:nbatch){
    DG   <- rbind(DG,
                  readRDS(file = paste0("./out/DG_",     var, "_", tag, "_", "batch=", batch, ".rds"))[,,m])
  }
  DGq[1, , m] <- apply(X = DG, MARGIN = 2, FUN = quantile, probs = 0.025)
  DGq[2, , m] <- apply(X = DG, MARGIN = 2, FUN = quantile, probs = 0.500)
  DGq[3, , m] <- apply(X = DG, MARGIN = 2, FUN = quantile, probs = 0.975)
  remove(DG)
  gc()
}

# Saves the Values
saveRDS(object = DAq,  file = paste0("./out/DA_q_", var, "_", tag, ".rds"))
saveRDS(object = DGq,  file = paste0("./out/DG_q_",  var, "_", tag, ".rds"))

# Saves g
# First Batch
g    <- readRDS(file = paste0("./out/g_",     var, "_", tag, "_", "batch=", "1", ".rds"))
# Rest of Batches
for(batch in 2:nbatch){
  g    <- rbind(g,
                readRDS(file = paste0("./out/g_",     var, "_", tag, "_", "batch=", batch, ".rds")))
}
g <- colMeans(g)
# Saves the Values
saveRDS(object = g,  file = paste0("./out/g_q_", var, "_", tag, ".rds"))

# Saves s2
# First Batch
s2    <- readRDS(file = paste0("./out/s2_",     var, "_", tag, "_", "batch=", "1", ".rds"))
# Rest of Batches
for(batch in 2:nbatch){
  s2    <- c(s2,
             readRDS(file = paste0("./out/s2_",     var, "_", tag, "_", "batch=", batch, ".rds")))
}
# Saves the Values
saveRDS(object = s2,  file = paste0("./out/s2_q_", var, "_", tag, ".rds"))
