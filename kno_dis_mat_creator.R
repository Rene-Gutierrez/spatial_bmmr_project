################################################################################
# Grid and Distance Matrix Computations
################################################################################
# SPN Region Selection
spn <- c(0, 1, 2, 3)
tag <- c()
for(i in spn){
  tag <- paste0(tag, i)
}

K <- 125

# Gets the Locations
L <- readRDS(file = paste0("./dat/rea_L_", tag, ".rds"))
C <- readRDS(file = paste0("./dat/rea_C_", tag, ".rds"))
P <- dim(C)[2]

max_R      <- matrix(data = NA, nrow = 3, ncol = P)
min_R      <- matrix(data = NA, nrow = 3, ncol = P)
max_R[1, ] <- apply(X = L[1,,], MARGIN = 2, FUN = max, na.rm = TRUE) 
min_R[1, ] <- apply(X = L[1,,], MARGIN = 2, FUN = min, na.rm = TRUE)
max_R[2, ] <- apply(X = L[2,,], MARGIN = 2, FUN = max, na.rm = TRUE) 
min_R[2, ] <- apply(X = L[2,,], MARGIN = 2, FUN = min, na.rm = TRUE)
max_R[3, ] <- apply(X = L[3,,], MARGIN = 2, FUN = max, na.rm = TRUE) 
min_R[3, ] <- apply(X = L[3,,], MARGIN = 2, FUN = min, na.rm = TRUE)

K3 <- round(K^(1/3))
kno <- list()
for(p in 1:P){
  print(p)
  knots_pos <- matrix(data = NA, nrow = K3, ncol = 3)
  for(i in 1:3){
    knots_pos[, i] <- seq(min_R[i, p], max_R[i, p], length.out = K3)
  }
  knots <- matrix(nrow = K, ncol = 3)
  k <- 1
  for(i in 1:K3){
    for(j in 1:K3){
      for(l in 1:K3){
        knots[k, 1] <- knots_pos[i, 1]
        knots[k, 2] <- knots_pos[j, 2]
        knots[k, 3] <- knots_pos[l, 3]
        k           <- k + 1
      }
    }
  }
  
  for(k in 1:K){
    while(sum((round(knots[k,])[1] == L[1,,p]) * (round(knots[k,])[2] == L[2,,p]) * (round(knots[k,])[3] == L[3,,p]), na.rm = TRUE) == 0){
      knots[k, 1] <- runif(n = 1, min = min_R[1, p], max = max_R[1, p])
      knots[k, 2] <- runif(n = 1, min = min_R[2, p], max = max_R[2, p])
      knots[k, 3] <- runif(n = 1, min = min_R[3, p], max = max_R[3, p])
    }
  }
  
  kno[[p]] <- knots
}

# Saves the Knots
saveRDS(object = kno, file = paste0("./dat/knots_", k, "_", tag, ".rds"))

# Creates the Distance Matrices
Dkk <- list()
Dks <- list()
for(p in 1:P){
  # Dkk
  D <- matrix(data = 0, nrow = K, ncol = K)
  for(k in 1:K){
    for(kk in 1:K){
      D[k, kk] <- sum((kno[[p]][k, ] - kno[[p]][kk, ])^2)
    }
  }
  D        <- sqrt(D)
  Dkk[[p]] <- D
  # Dks
  loc <- t(L[,C[, p], p])
  N   <- dim(loc)[1]
  D <- matrix(data = 0, nrow = K, ncol = N)
  for(k in 1:K){
    for(n in 1:N){
      D[k, n] <- sum((kno[[p]][k, ] - loc[n, ])^2)
    }
  }
  D   <- sqrt(D)
  Dks[[p]] <- D  
}
# Saves the Distance Matrices
saveRDS(object = Dkk, file = paste0("./dat/Dkk_", k, "_", tag, ".rds"))
saveRDS(object = Dks, file = paste0("./dat/Dks_", k, "_", tag, ".rds"))

for(p in 1:P){
  print(max(apply(X = Dks[[p]], MARGIN = 2, FUN = min)))
}
