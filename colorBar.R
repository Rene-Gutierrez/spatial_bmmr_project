colorBar <- function(hor = TRUE,
                     txt = "",
                     vmax,
                     vmin = -vmax,
                     cex  = 1){
  if(hor){
    par(bg = rgb(1, 1, 1))
    plot(NULL,
         ylim = c(0, 1),
         xlim = c(-1, 1),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         yaxt = 'n',
         xaxs = 'i',
         yaxs = 'i')
    sca    <- seq(-1, 1, length.ou = 101)
    scaDif <- sca[2] - sca[1] 
    rect(xleft   = sca,
         ybottom = 0,
         xright  = sca + scaDif,
         ytop    = 1,
         col     = rgb((1 + sign(sca)) / 2, 0, abs(sign(sca) - 1) / 2, abs(sca)),
         border  = rgb((1 + sign(sca)) / 2, 0, abs(sign(sca) - 1) / 2, abs(sca)))
    if(vmin != -vmax){
      axis(side = 1, at = c(-1, 0, 1), labels = round(c( vmin, 0, vmax), 2))
    } else {
      axis(side = 1, at = c(-1, 0, 1), labels = round(c(-vmax, 0, vmax), 2))
    }
    mtext(text = txt, side = 1, line = 3, cex = cex)
  } else {
    par(bg = rgb(1, 1, 1))
    plot(NULL,
         ylim = c(-1, 1),
         xlim = c(0, 1),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         yaxt = 'n',
         xaxs = 'i',
         yaxs = 'i')
    sca    <- seq(-1, 1, length.ou = 101)
    scaDif <- sca[2] - sca[1]
    rect(xleft   = 0,
         ybottom = sca,
         xright  = 1,
         ytop    = sca + scaDif,
         col     = rgb((1 + sign(sca)) / 2, 0, abs(sign(sca) - 1) / 2, abs(sca)),
         border  = rgb((1 + sign(sca)) / 2, 0, abs(sign(sca) - 1) / 2, abs(sca)))
    if(vmin != -vmax){
      axis(side = 4, at = c(-1, 0, 1), labels = round(c( vmin, 0, vmax), 2))
    } else {
      axis(side = 4, at = c(-1, 0, 1), labels = round(c(-vmax, 0, vmax), 2))
    }
    
    mtext(text = txt, side = 4, line = 2, las = 2, cex = cex)
  }
}