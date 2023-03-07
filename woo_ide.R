################################################################################
# Woodbury Matrix Identity
################################################################################

# Inputs
# iA: Inverse Diagonal Matrix
# U:  Side Matrices
# iC: Inverse Regular Matrix

woo_ide <- function(iA,
                    U,
                    iC){
  # Auxiliary Computations
  M   <- iC + U %*% (diag(iA) * t(U))
  AU  <- diag(iA) * t(U)
  UMU <- AU %*% solve(M, t(AU))
  
  # Returns the Inverse
  return(iA - UMU)
}