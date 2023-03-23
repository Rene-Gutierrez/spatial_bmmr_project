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

# Data Import
G      <- readRDS(file = paste0("./dat/rea_G_", tag, ".rds"))
C      <- readRDS(file = paste0("./dat/rea_C_", tag, ".rds"))
L      <- readRDS(file = paste0("./dat/rea_L_", tag, ".rds"))
reg_id <- readRDS(file = paste0("./dat/regID_", tag, ".rds"))
mG     <- apply(X = G, MARGIN = c(2, 3), FUN = mean)

# List
P  <- dim(L)[3]
V  <- dim(L)[2]
lL <- list()
for(p in 1:P){
  lL[[p]] <- list()
  for(v in 1:V){
    lL[[p]][[v]] <- c(L[1, v, p], L[2, v, p], L[3, v, p])
  }
}

M  <- rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 1, 0), c(0, -1, 0), c(0, 0, 1), c(0, 0, -1))
LN <- list()
for(i in 1:6){
  LN[[i]] <- list()
  for(p in (1:10)[-c(3,4,5,6)]){
    LN[[i]][[p]] <- list()
    for(v in 1:sum(C[,p])){
      LN[[i]][[p]][[v]] <- c(L[1, v, p], L[2, v, p], L[3, v, p]) + M[i,]
    }
  }
}


GN <- array(data = NA, dim = c(6, V, P))
for(i in 1:6){
  print(i)
  for(p in (1:10)[-c(3,4,5,6)]){
    print(p)
    for(v in 1:sum(C[,p])){
      temp <- lL[[p]] %in% list(lL[[p]][[v]] + M[i,])
      if(sum(temp)){
        GN[i, v, p] <- mG[which(temp), p]
      }
    }
  }
}
