spa_iboom_iterator <- function(G,
                               A,
                               y,
                               X,
                               XG,
                               Zv,
                               Theta,
                               B,
                               DA,
                               DG,
                               wG,
                               C,
                               t2T,
                               l2T,
                               vT,
                               xiT,
                               t2B,
                               l2B,
                               vB,
                               xiB,
                               s2,
                               g,
                               Dkk,
                               Dks,
                               ga = 0.5 * s2,
                               ph = 1){
  # Problem dimensions
  N  <- length(y)            # Number of Observations
  mV <- dim(B)[1]            # Maximum Voxel Size
  V  <- colSums(C)           # Number of Voxels per Region
  P  <- dim(B)[2]            # Number of ROI's
  Q  <- sum(g)               # Number of Active Regions
  M  <- ncol(X)              # Number of Covariates
  gp <- numeric(length = P)
  
  # Samples DA
  for(p in 1:P){
    Z <- matrix(data = NA, nrow = 0, ncol = M)
    R <- c(A[,-p, p]) - kronecker(X = Theta[-p, p], y)
    for(pp in (1:P)[-p]){
      Z <- rbind(Z, t(t(X) * DA[pp, ]))
    }
    ZZ      <- solve(t(Z) %*% Z)
    DAphat  <- ZZ %*% t(Z) %*% R
    DA[p, ] <- c(DAphat) + rmvnorm(n = 1, sigma = s2 * ZZ)
  }
  # Samples DG
  for(p in 1:P){
    R       <- c(t(t(G[,,p]) - wG[, p])) - kronecker(X = B[, p], Y = y)
    Z       <- XG[!is.na(R), ]
    R       <- R[!is.na(R)]
    ZZ      <- solve(t(Z) %*% Z)
    DGphat  <- ZZ %*% t(Z) %*% R
    DG[p, ] <- c(DGphat) + rmvnorm(n = 1, sigma = s2 * ZZ)
  }
  # Creates DT and DB
  DT   <- array(data = 0,  dim = c(M, P, P))
  DB   <- array(data = NA, dim = c(M, mV, P))
  D    <- matrix(data = NA, nrow = M, ncol = (P * (P - 1) / 2 + sum(C)))
  for(m in 1:M){
    DT[m, , ]       <- DA[, m] %*% t(DA[, m])
    diag(DT[m, , ]) <- NA
    DB[m, , ][C]    <- kronecker(X = t(DG[, m]), Y = rep(1, mV))[C]
    D[m, ]          <- c(DT[m,,][upper.tri(DT[m,,])], DB[m,,][C])
  }
  
  # Samples w_tilde
  AP <- A - apply(DT, c(2, 3), function(x) X %*% x)
  GP <- G - apply(DB, c(2, 3), function(x) X %*% x)
  wG <- matrix(data = NA, nrow = mV, ncol = P)
  for(p in 1:P){
    yb <- colMeans(GP[,,p] - kronecker(X = y, Y = t(B[, p])))[C[,p]]
    wt <- w_tilde_sampler(s2  = s2,
                          ga  = ga,
                          ph  = ph,
                          Dkk = Dkk[[p]],
                          Dks = Dks[[p]],
                          yb  = yb,
                          K   = N)
    wG[C[,p],p] <- wt
  }
  
  # Samples g, Theta and B
  GP  <- GP - array(data = wG %x% rep(1, N), dim = c(N, mV, P))
  Say <- apply(X = (AP * y), MARGIN = c(2, 3), sum)
  Sgy <- apply(X = (GP * y), MARGIN = c(2, 3), sum)
  Syy <- sum(y^2)
  for(p in 1:P){
    res <- group_iterator(Say = Say, 
                          Sgy = Sgy,
                          Syy = Syy,
                          C   = C,
                          LT  = t2T * l2T,
                          LB  = t2B * l2B,
                          s2  = s2,
                          g   = g,
                          p   = p)
    # Updates gp
    gp[p] <- res$pr
    # Updates g
    g <- res$g
    # Updates Theta
    Q <- sum(g[-p])
    if(Q > 0){
      Theta[-p, p][g[-p] == 1] <- res$b[1:Q]
      Theta[p, -p][g[-p] == 1] <- Theta[-p, p][g[-p] == 1]
    }
    # Updates B
    if(Q > 0){
      B[C[, p], p] <- res$b[-(1:Q)]
    } else {
      B[C[, p], p] <- res$b
    }
  }
  
  # Horseshoe Structure for B
  # Samples l2B
  for(p in 1:P){
    if(g[p] == 1){
      l2B[1:V[p], p] <- 1 / rgamma(n     = V[p],
                                   shape = 1,
                                   rate  = 1 / vB[, p] + B[, p]^2 / (2 * s2 * t2B))
    } else {
      l2B[1:V[p], p] <- 1 / rgamma(n     = V[p],
                                   shape = 1 / 2,
                                   rate  = 1 / vB[, p])
    }
  }
  
  # Samples t2B
  if(sum(g) > 0){
    t2B <- 1 / rgamma(n     = 1,
                      shape = (sum(V[g == 1]) + 1) / 2,
                      rate  = 1 / xiB +
                        sum(B^2 / (2 * s2 * l2B), na.rm = TRUE))
  } else {
    t2B <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiB)
  }
  
  # Samples vB
  for(p in 1:P){
    vB[1:V[p], p]<- 1 / rgamma(n     = V[p],
                               shape = 1,
                               rate  = 1 + 1 / l2B[V[p], p])
  }
  
  # Samples xiB
  xiB <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2B)
  
  # Horseshoe Structure for Theta
  # Samples l2T
  for(p in 2:P){
    for(pp in 1:(p - 1)){
      if(g[p] * g[pp] == 1){
        l2T[p, pp] <- 1 / rgamma(n     = 1,
                                 shape = 1,
                                 rate  = 1 / vT[p, pp] + Theta[p, pp]^2 / (2 * s2 * t2T))
        l2T[pp, p] <- l2T[p, pp]
      } else {
        l2T[p, pp] <- 1 / rgamma(n     = 1,
                                 shape = 1 / 2,
                                 rate  = 1 / vT[p, pp])
        l2T[pp, p] <- l2T[p, pp]
      }
    }
  }
  
  # Samples t2T
  if(sum(g) > 0){
    t2T <- 1 / rgamma(n     = 1,
                      shape = (sum(g) * (sum(g) - 1) / 2 + 1) / 2,
                      rate  = 1 / xiT +
                        sum(Theta[lower.tri(Theta)]^2 / (2 * s2 * l2T[lower.tri(l2T)])))
  } else {
    t2T <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiT)
  }
  
  # Samples vT
  out <- 1 / rgamma(n     = P * (P - 1) / 2,
                    shape = 1,
                    rate  = 1 + 1 / l2T[lower.tri(l2T)])
  vT[upper.tri(vT, diag = TRUE)] <- 0
  vT[lower.tri(vT)]              <- out
  vT                             <- vT + t(vT)
  
  # Samples xiT
  xiT <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2T)
  
  # Samples s2
  Q    <- sum(g == 1)
  Beta <- c(Theta[upper.tri(Theta)], B[C])
  bs2  <- 2   + sum((Zv - X %*% D - y %*% t(Beta))^2)
  bs2  <- bs2 + sum(Theta[g == 1, g == 1]^2 / l2T[g == 1, g == 1], na.rm = TRUE) / 2 / t2T
  bs2  <- bs2 + sum(B[, g == 1]^2 / l2B[, g == 1], na.rm = NA) / t2B
  bs2  <- bs2 / 2
  as2  <- 2 + N * P * (P - 1) / 2
  as2  <- as2 + N * sum(C)
  as2  <- as2 + Q * (Q - 1) / 2
  as2  <- as2 + sum(C[, g == 1])
  as2  <- as2 / 2
  s2   <- 1 / rgamma(n = 1, shape = as2, rate = bs2)
  
  # Returns Values
  return(list(Theta = Theta,
              B     = B,
              DA    = DA,
              DG    = DG,
              wG    = wG,
              g     = g,
              s2    = s2,
              as2   = as2,
              bs2   = bs2,
              l2T   = l2T,
              t2T   = t2T,
              vT    = vT,
              xiT   = xiT,
              l2B   = l2B,
              t2B   = t2B,
              vB    = vB,
              xiB   = xiB,
              gp    = gp,
              wG    = wG))
}