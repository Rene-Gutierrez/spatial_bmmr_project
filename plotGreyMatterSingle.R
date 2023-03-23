plotGreyMatterSingle <- function(L,
                                 C,
                                 D,
                                 z,
                                 reg1,
                                 reg2,
                                 line_w,
                                 title = NULL,
                                 vmax){
  # Plotting Dimensions
  xmax <- max(L[1,,][C])
  xmin <- min(L[1,,][C])
  ymax <- max(L[2,,][C])
  ymin <- min(L[2,,][C])
  
  # Plotting Area
  if(is.null(title)){
    plot(NULL,
         type = 'l',
         ylim = c(ymin, ymax),
         xlim = c(xmin, xmax),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         yaxt = 'n')
  } else {
    plot(NULL,
         type = 'l',
         ylim = c(ymin, ymax),
         xlim = c(xmin, xmax),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         yaxt = 'n',
         main = title)
  }
  rect(xleft   = par("usr")[1],
       ybottom = par("usr")[3],
       xright  = par("usr")[2],
       ytop    = par("usr")[4],
       col     = rgb(0, 0, 0, 0.2),
       border  = rgb(0, 0, 0, 0.2))
  text(x = xmin + 10,
       y = ymin + 1,
       labels = c(paste0("z = ", z)))
  
  # Loops Trough every Region
  for(p in reg1){
    x   <- L[1, C[, p], p][L[3, C[, p], p] == z]
    y   <- L[2, C[, p], p][L[3, C[, p], p] == z]
    val <- D[C[, p], p][L[3, C[, p], p] == z]
    
    if(length(val) > 0){
      ux <- unique(x)
      uy <- unique(y)
      nx <- length(ux)
      ny <- length(uy)
      
      # Region Boundaries
      lx_min <- list()
      lx_max <- list()
      j <- 1
      for(i in ux){
        lx_min[[j]] <- y[x == i][!(y[x == i] %in% (y[x == i] + 1))]
        lx_max[[j]] <- y[x == i][!(y[x == i] %in% (y[x == i] - 1))]
        j <- j + 1
      }
      
      ly_min <- list()
      ly_max <- list()
      j <- 1
      for(i in uy){
        ly_min[[j]] <- x[y == i][!(x[y == i] %in% (x[y == i] + 1))]
        ly_max[[j]] <- x[y == i][!(x[y == i] %in% (x[y == i] - 1))]
        j <- j + 1
      }
      
      for(i in 1:nx){
        segments(x0  = ux[i],
                 y0  = lx_min[[i]],
                 x1  = ux[i] + 1,
                 y1  = lx_min[[i]],
                 lwd = line_w)
        segments(x0  = ux[i],
                 y0  = lx_max[[i]] + 1,
                 x1  = ux[i] + 1,
                 y1  = lx_max[[i]] + 1,
                 lwd = line_w)
      }
      
      for(i in 1:ny){
        segments(x0  = ly_min[[i]],
                 y0  = uy[i],
                 x1  = ly_min[[i]],
                 y1  = uy[i] + 1,
                 lwd = line_w)
        segments(x0  = ly_max[[i]] + 1,
                 y0  = uy[i],
                 x1  = ly_max[[i]] + 1,
                 y1  = uy[i] + 1,
                 lwd = line_w)
      }
      
      if(p %in% reg2){
        rect(xleft   = x,
             ybottom = y,
             xright  = x + 1,
             ytop    = y + 1,
             col     = rgb(1, 1, 1),
             border  = rgb(1, 1, 1))
        # rect(xleft   = x,
        #      ybottom = y,
        #      xright  = x + 1,
        #      ytop    = y + 1,
        #      col     = rgb((1 + sign(val)) / 2, 0, abs(sign(val) - 1) / 2, abs(val)),
        #      border  = rgb((1 + sign(val)) / 2, 0, abs(sign(val) - 1) / 2, abs(val)))
        rect(xleft   = x,
             ybottom = y,
             xright  = x + 1,
             ytop    = y + 1,
             col     = rgb((1 + sign(val)) / 2, 0, abs(sign(val) - 1) / 2, abs(val)),
             border  = NA)
      } else {
        # rect(xleft   = x,
        #      ybottom = y,
        #      xright  = x + 1,
        #      ytop    = y + 1,
        #      col     = rgb(0, 0, 0, 0.10),
        #      border  = rgb(0, 0, 0, 0.10))
        rect(xleft   = x,
             ybottom = y,
             xright  = x + 1,
             ytop    = y + 1,
             col     = rgb(1, 1, 1),
             border  = rgb(1, 1, 1))
        rect(xleft   = x,
             ybottom = y,
             xright  = x + 1,
             ytop    = y + 1,
             col     = rgb(0, 0, 0, 0.10),
             border  = NA)
      }
    }
  }
}