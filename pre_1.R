################################################################################
# Real Data Analysis
################################################################################

### Seed
set.seed(25052022)
### Outcome
var <- 48

### SPN Selection
spn <- c(0, 1, 2, 3)
tag <- c()
for(i in spn){
  tag <- paste0(tag, i)
}
reg      <- readRDS(file = paste0("./dat/regID_", tag, ".rds"))
reg_full <- reg

### Data Import
A      <- readRDS(file = paste0("./dat/rea_A_", tag, ".rds"))
G      <- readRDS(file = paste0("./dat/rea_G_", tag, ".rds"))
y      <- readRDS(file =        "./dat/rea_y.rds")
C      <- readRDS(file = paste0("./dat/rea_C_", tag, ".rds"))
L      <- readRDS(file = paste0("./dat/rea_L_", tag, ".rds"))
reg_id <- readRDS(file = paste0("./dat/regID_", tag, ".rds"))

### Data Wrangling and Variable Selection
var <- 48
X   <- y[, c(5, 12, 13)]
y   <- y[, var]
# Removes NA
selX <- rowSums(is.nan(as.matrix(X))) == 0
sel  <- !is.na(y)
sel  <- sel & selX
y    <- y[sel]
X    <- X[sel,]
A    <- A[sel, ,]
G    <- G[sel, ,]

### Dimensions
N  <- length(y)
P  <- dim(G)[3]
mV <- dim(G)[2]

### Normalizes
# Centers
mA <- apply(X = A, MARGIN = c(2, 3), FUN = mean)
mG <- apply(X = G, MARGIN = c(2, 3), FUN = mean)
A  <- A - array(data = mA %x% rep(1, N), dim = c(N, P,  P))
# G  <- G - array(data = mG %x% rep(1, N), dim = c(N, mV, P))
y  <- (y - mean(y)) / sd(y)

Zv <- matrix(data = NA, nrow = N, ncol = sum(C) + P * (P - 1) / 2)
for(i in 1:N){
  Zv[i, ] <- c(A[i,,][upper.tri(A[i,,])], G[i,,][C])
}

# Generates Covariates
X[, 1] <- (X[, 1] == 1) + 0
X[, 2] <- X[, 2] - mean(X[, 2])
X[, 3] <- X[, 3] - mean(X[, 3])
X      <- as.matrix(X)

XG  <- kronecker(X = rep(1, mV), Y = X)

saveRDS(object = X, file = paste0("./dat/rea_X_", var, ".rds"))
saveRDS(object = y, file = paste0("./dat/rea_y_", var, ".rds"))
saveRDS(object = A, file = paste0("./dat/rea_A_", var, "_", tag, ".rds"))
saveRDS(object = G, file = paste0("./dat/rea_G_", var, "_", tag, ".rds"))