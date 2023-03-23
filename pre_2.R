### Libraries
source("./colorBar.R")
source("./plotGreyMatterSingle.R")

# Patient
reg <- (1:10)[-c(3,4,5,6)]

# Patient
pat <- 1
mG  <- apply(X = GN, MARGIN = c(2, 3), FUN = mean, na.rm = TRUE)
mG  <- mG - apply(X = G, MARGIN = c(2, 3), FUN = mean)
eG  <- matrix(data = NA, nrow = dim(mG)[1], ncol = dim(mG)[2])
eG[C] <- rnorm(n = sum(C), sd = 0.05)
mG  <- mG + eG

# Z level
uni <- unique(L[3,,][C])
sli <- 8
zs  <- round(seq(min(uni), max(uni), length.out = sli + 2))[-c(1, sli + 2)]

# Difference Between Covariates and Covariate-less
D1   <- G[pat, ,] - apply(X = G, MARGIN = c(2, 3), FUN = mean)
D1[, -reg] <- NA
# D1   <- mG
vmax <- max(D1, na.rm = TRUE)
vmin <- min(D1, na.rm = TRUE)

# Maximum Values
D1 <- (D1 - vmin) / (vmax - vmin)
D1 <- D1 * 2 - 1

# Open a pdf file
filNam <- paste0("./plo/cov_2_grey_matter_patient_", pat, "_", var, "_", tag, ".pdf")
pdf(file   = filNam,
    width  = 6,
    height = 4.2)
# 2. Create a plot
#par(bg = rgb(0, 0, 0, 0.25))
par(mar = c(0.2, 0.2, 0.2, 0.2))
layout(mat = matrix(data  = c(1, 2, 3, 4,
                              5, 6, 7, 8,
                              9, 9, 9, 9),
                    nrow  = 3,
                    ncol  = 4,
                    byrow = TRUE),
       heights = c(10, 10, 4),
       widths  = c(1, 1, 1, 1))
for(k in 1:sli){
  if(k / sli <= 0.5){
    par(mar = c(0.2, 0.2, 2, 0.2))
  } else {
    par(mar = c(0.2, 0.2, 0.2, 0.2))
  }
  z <- zs[k]
  if(k == 1){
    plotGreyMatterSingle(L    = L,
                         C    = C,
                         D    = D1,
                         z    = z,
                         reg1 = reg_full,
                         reg2 = reg,
                         line_w = 3,
                         title = "")
  } else {
    plotGreyMatterSingle(L    = L,
                         C    = C,
                         D    = D1,
                         z    = z,
                         reg1 = reg_full,
                         reg2 = reg,
                         line_w = 3)
  }
}
# Color Bar
par(bg = rgb(1, 1, 1))
par(mar = c(4, 1, 0, 1))
colorBar(vmax = vmax,
         vmin = vmin,
         txt  = "w")
# Close the pdf file
dev.off()

# Z level
uni <- unique(L[3,,][C])
sli <- 8
zs  <- round(seq(min(uni), max(uni), length.out = sli + 2))[-c(1, sli + 2)]

# Difference Between Covariates and Covariate-less
# D1   <- G[pat, ,]
# D1   <- apply(X = G, MARGIN = c(2, 3), FUN = mean)
# D1   <- mG
vmax <- max(D1, na.rm = TRUE)
vmin <- min(D1, na.rm = TRUE)

# Maximum Values
D1 <- (D1 - vmin) / (vmax - vmin)
D1 <- D1 * 2 - 1

for(k in 1:sli){
  z <- zs[k]
  # Open a pdf file
  filNam <- paste0("./plo/cov_pro_grey_matter_patient_", pat, "_", z, "_", var, "_", tag, ".pdf")
  pdf(file   = filNam,
      width  = 4,
      height = 6.2)
  # 2. Create a plot
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  layout(mat = matrix(data  = c(1, 2),
                      nrow  = 2,
                      ncol  = 1,
                      byrow = TRUE),
         heights = c(5, 1),
         widths  = c(1, 1))
  plotGreyMatterSingle(L      = L,
                       C      = C,
                       D      = D1,
                       z      = z,
                       reg1   = reg_full,
                       reg2   = reg,
                       line_w = 3,
                       title  = "")
  # Color Bar
  par(bg = rgb(1, 1, 1))
  par(mar = c(4, 0.2, 0.2, 0.2))
  colorBar(vmax = vmax,
           vmin = vmin,
           txt  = "w",
           cex  = 1.5)
  # Close the pdf file
  dev.off()
}