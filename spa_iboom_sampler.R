spa_iboom_sampler <- function(y,
                              X,
                              XG,
                              Zv,
                              G,
                              A,
                              C,
                              Dkk,
                              Dks,
                              nmcmc  = 1000,
                              burnin = 0,
                              onlyg  = TRUE,
                              state  = list()){
  # Computes Sufficient Statistics and Problem Dimensions
  N   <- dim(G)[1]
  mV  <- dim(G)[2]
  P   <- dim(G)[3]
  M   <- dim(X)[2]
  
  # Sample Holders
  if(onlyg){
    sg     <- matrix(data    = NA, nrow = nmcmc, ncol = P)
    sB     <- array(data     = NA, dim  = c(nmcmc, P * (P - 1)/ 2 + sum(C)))
    swG    <- array(data     = NA, dim  = c(nmcmc, mV, P))
    sDA    <- array(data     = NA, dim  = c(nmcmc, P, M))
    sDG    <- array(data     = NA, dim  = c(nmcmc, P, M))
    ss2    <- numeric(length = nmcmc)
  } else {
    sTheta <- array(data     = NA, dim  = c(burnin + nmcmc, P, P))
    sB     <- array
    sDA    <- array(data     = NA, dim  = c(nmcmc, P, M))
    sDG    <- array(data     = NA, dim  = c(nmcmc, P, M))
    sgp    <- matrix(data    = NA, nrow = burnin + nmcmc, ncol = P)
    ss2    <- numeric(length = burnin + nmcmc)
    sl2T   <- array(data     = NA, dim  = c(burnin + nmcmc, P, P))
    st2T   <- numeric(length = burnin + nmcmc)
    svT    <- array(data     = NA, dim  = c(burnin + nmcmc, P, P))
    sxiT   <- numeric(length = burnin + nmcmc)
    sl2B   <- array(data     = NA, dim  = c(burnin + nmcmc, mV, P))
    st2B   <- numeric(length = burnin + nmcmc)
    svB    <- array(data     = NA, dim  = c(burnin + nmcmc, mV, P))
    sxiB   <- numeric(length = burnin + nmcmc)
  }
  
  # Initialization
  if(length(state) != 0){
    Theta  <- state$Theta
    B      <- state$B
    wG     <- state$wG
    DA     <- state$DA
    DG     <- state$DG
    t2T    <- state$t2T
    l2T    <- state$l2T
    xiT    <- state$xiT 
    vT     <- state$vT
    t2B    <- state$t2B
    l2B    <- state$l2B
    xiB    <- state$xiB
    vB     <- state$vB
    s2     <- state$s2
    g      <- state$g
  } else {
    Theta  <- matrix(data = 0,  nrow = P,  ncol = P)
    B      <- matrix(data = NA, nrow = mV, ncol = P)
    B[C]   <- 0
    wG     <- matrix(data = NA, nrow = mV, ncol = P)
    wG[C]  <- 0
    DA     <- matrix(data = rnorm(n = P * M,
                                  mean = 0,
                                  sd   = 1/10),
                     nrow = P,  ncol = M)
    DG     <- matrix(data = 0, nrow = P,  ncol = M)
    t2T    <- 1
    l2T    <- matrix(data = 1, nrow = P,   ncol = P)
    xiT    <- 1
    vT     <- matrix(data = 1, nrow = P,   ncol = P)
    t2B    <- 1
    l2B    <- matrix(data = NA, nrow = mV, ncol = P)
    l2B[C] <- 1
    xiB    <- 1
    vB     <- matrix(data = NA, nrow = mV, ncol = P)
    vB[C]  <- 1
    s2     <- 1
    g      <- rep(1, P)
  }
  
  # Progress Bar
  pb <- txtProgressBar(min     = 0,
                       max     = 1,
                       initial = 0,
                       style   = 3,
                       width   = 72)
  
  # First Run
  out <- spa_iboom_iterator(G     = G,
                            A     = A,
                            y     = y,
                            X     = X,
                            XG    = XG,
                            Zv    = Zv,
                            Theta = Theta,
                            B     = B,
                            DA    = DA,
                            DG    = DG,
                            wG    = wG,
                            C     = C,
                            t2T   = t2T,
                            l2T   = l2T,
                            vT    = vT,
                            xiT   = xiT,
                            t2B   = t2B,
                            l2B   = l2B,
                            vB    = vB,
                            xiB   = xiB,
                            s2    = s2,
                            g     = g,
                            Dkk   = Dkk,
                            Dks   = Dks)
  
  # Sampling
  for(s in 1:(burnin + nmcmc)){
    # Runs Iterator
    out <- spa_iboom_iterator(G     = G,
                              A     = A,
                              y     = y,
                              X     = X,
                              XG    = XG,
                              Zv    = Zv,
                              Theta = out$Theta,
                              B     = out$B,
                              DA    = out$DA,
                              DG    = out$DG,
                              wG    = out$wG,
                              C     = C,
                              t2T   = out$t2T,
                              l2T   = out$l2T,
                              vT    = out$vT,
                              xiT   = out$xiT,
                              t2B   = out$t2B,
                              l2B   = out$l2B,
                              vB    = out$vB,
                              xiB   = out$xiB,
                              s2    = out$s2,
                              g     = out$g,
                              Dkk   = Dkk,
                              Dks   = Dks)
    # Saves Samples
    if(s > burnin){
      if(onlyg){
        sg[s - burnin,]   <- out$g
        sB[s - burnin,]   <- c(out$Theta[upper.tri(out$Theta)], out$B[C])
        swG[s - burnin,,] <- out$wG
        sDA[s - burnin,,] <- out$DA
        sDG[s - burnin,,] <- out$DG
        ss2[s - burnin]   <- out$s2
      }else {
        sTheta[s,,] <- out$Theta
        sB[s,,]     <- out$B
        sgp[s,]     <- out$gp
        ss2[s]      <- out$s2
        sl2T[s,,]   <- out$l2T
        st2T[s]     <- out$t2T
        svT[s,,]    <- out$vT
        sxiT[s]     <- out$xiT
        sl2B[s,,]   <- out$l2B
        st2B[s]     <- out$t2B
        svB[s,,]    <- out$vB
        sxiB[s]     <- out$xiB
      }
    }
    
    # Progress Bar Update
    setTxtProgressBar(pb    = pb,
                      value = s / (burnin + nmcmc))
  }
  
  # Splits Burn-In and After-Burin Samples
  if(onlyg){
    # sam <- list(g     = sg[(burnin + 1):(burnin + nmcmc),],
    #             Theta = sTheta[(burnin + 1):(burnin + nmcmc),,],
    #             B     = sB[(burnin + 1):(burnin + nmcmc),,],
    #             s2    = ss2[(burnin + 1):(burnin + nmcmc)])
    sam <- list(g     = sg,
                B     = sB,
                wG    = swG,
                DA    = sDA,
                DG    = sDG,
                s2    = ss2)
  } else {
    sam <- list(Theta = sTheta[(burnin + 1):(burnin + nmcmc),,],
                B     = sB[(burnin + 1):(burnin + nmcmc),,],
                DA    = sDA[(burnin + 1):(burnin + nmcmc),,],
                DG    = sDG[(burnin + 1):(burnin + nmcmc),,],
                t2T   = st2T[(burnin + 1):(burnin + nmcmc)],
                l2T   = sl2T[(burnin + 1):(burnin + nmcmc),,],
                xiT   = sxiT[(burnin + 1):(burnin + nmcmc)],
                vT    = svT[(burnin + 1):(burnin + nmcmc),,],
                t2B   = st2B[(burnin + 1):(burnin + nmcmc)],
                l2B   = sl2B[(burnin + 1):(burnin + nmcmc),,],
                xiB   = sxiB[(burnin + 1):(burnin + nmcmc)],
                vB    = svB[(burnin + 1):(burnin + nmcmc),,],
                s2    = ss2[(burnin + 1):(burnin + nmcmc)],
                g     = sg[(burnin + 1):(burnin + nmcmc),],
                gp    = sgp[(burnin + 1):(burnin + nmcmc),])
  }
  
  # Saves the State
  state = list(Theta = out$Theta,
               B     = out$B,
               wG    = out$wG,
               DA    = out$DA,
               DG    = out$DG,
               t2T   = out$t2T,
               l2T   = out$l2T,
               xiT   = out$xiT,
               vT    = out$vT,
               t2B   = out$t2B,
               l2B   = out$l2B,
               xiB   = out$xiB,
               vB    = out$vB,
               s2    = out$s2,
               g     = out$g)
  
  # Returns Values
  return(list(sam   = sam,
              state = state,
              nmcmc = nmcmc))
}