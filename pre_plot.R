### Libraries
source("./colorBar.R")
source("./plotGreyMatterSingle.R")

# Region Selection
# Full Region
spn <- c(0, 1, 2, 3)
tag <- c()
for(i in spn){
  tag <- paste0(tag, i)
}
C       <- readRDS(file = paste0("./dat/rea_C_", tag, ".rds"))
L       <- readRDS(file = paste0("./dat/rea_L_", tag, ".rds"))
P       <- dim(C)[2]
mV      <- dim(C)[1]
reg_full <- readRDS(file = paste0("./dat/regID_", tag, ".rds"))
# Full Region
spn <- c(0, 1, 2, 3)
tag <- c()
for(i in spn){
  tag <- paste0(tag, i)
}
reg      <- readRDS(file = paste0("./dat/regID_", tag, ".rds"))

# Recovers the Statistics
# Beta    <- readRDS(file = paste0("./out/Beta_q_",  var, "_", tag, ".rds"))
mV1     <- dim(C)[1]
# g       <- readRDS(file = paste0("./out/g_q_",     var, "_", tag, ".rds"))
g <- colMeans(out_rea$sam$g)

# Z level
uni <- unique(L[3,,][C])
sli <- 8
zs  <- round(seq(min(uni), max(uni), length.out = sli + 2))[-c(1, sli + 2)]

# Z level
uni <- unique(L[3,,][C])
sli <- 8
zs  <- round(seq(min(uni), max(uni), length.out = sli + 2))[-c(1, sli + 2)]

# Difference Between Covariates and Covariate-less
# D1 <- apply(X = out_rea$sam$wG, MARGIN = c(2, 3), FUN = mean, na.rm = TRUE)
# D1 <- matrix(data = NA, nrow = mV, ncol = P)
# D1[-((mV1 + 1):(mV)), reg] <- apply(X = out_rea$sam$wG, MARGIN = c(2, 3), FUN = mean, na.rm = TRUE)
D1   <- apply(X = G, MARGIN = c(2, 3), FUN = mean, na.rm = TRUE)
vmax <- max(abs(D1), na.rm = TRUE)
# vmin <- min(D1, na.rm = TRUE)
vmin <- -vmax

# Maximum Values
D1 <- (D1 - vmin) / (vmax - vmin)
D1 <- D1 * 2 - 1

# Open a pdf file
filNam <- paste0("./plo/spa_0_grey_matter_patient_", "_", var, "_", tag, ".pdf")
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
         txt  = "B")
# Close the pdf file
dev.off()

# Z level
uni <- unique(L[3,,][C])
sli <- 8
zs  <- round(seq(min(uni), max(uni), length.out = sli + 2))[-c(1, sli + 2)]

# Difference Between Covariates and Covariate-less
vmax <- max(D1, na.rm = TRUE)
vmin <- min(D1, na.rm = TRUE)

# Maximum Values
D1 <- (D1 - vmin) / (vmax - vmin)
D1 <- D1 * 2 - 1

for(k in 1:sli){
  z <- zs[k]
  # Open a pdf file
  filNam <- paste0("./plo/spa_0_grey_matter_patient_", "_", z, "_", var, "_", tag, ".pdf")
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