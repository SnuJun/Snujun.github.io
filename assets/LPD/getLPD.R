## input : (1) A: symmetric matrix and mostly not positive definite
##         (2) lwr_bd: a lower bound of eigenvalues of the output matrix
##         (3) type_norm: a type of matrix norms; 
##                        "spectral", "Frobenius", 
##                        "Linf", (column maximum norm), 
##                        "elemax" (elementwise maximum norm)
##         (4) mu: automatically chosen, unless specified. Default is NULL. 
##                 Do not specify it unless knowing what's going on.
##
## output: (1) mu, lower bound
##         (2) alpha, 
##         (3) LPD matrix (= alpha * A + mu * (1 - alpha) * A)
getLPD <- function(A, lwr_bd, type_norm, mu = NULL){
  
  res_EV <- getExtremeEV(A)
  if(is.null(mu)){
    switch(as.character(type_norm),
           "spectral" = {
             mu <- max(0.5 * (res_EV$largest + res_EV$smallest), 
                       lwr_bd)
             mu_range <- c(mu, Inf)
           },
           "Frobenius" = {
             res_eig <- eigen(A, only.values = TRUE)
             mu <- sum((res_eig$values - mean(res_eig$values))^2) / 
               sum((res_eig$values - min(res_eig$values))^2)
             mu_range <- c(mu, mu)
           },
           "Linf" = {
             M1 <- max(rowSums(abs(A)))
             A_tmp <- A - 2 * diag(diag(A))
             A_tmp[row(A_tmp) != col(A_tmp)] <- abs(A_tmp[row(A_tmp) != col(A_tmp)])
             M2 <- max(rowSums(A_tmp))
             
             if(res_EV$smallest + M2 > 1e-10){
               mu <- 1e+8
               mu_range <- c(Inf, Inf)
               warning(sprintf("The loss function is decreasing in mu.\nA very large mu %.1e is returned.", mu))
             } else if(res_EV$smallest + M2 < -1e-10){
               mu <- 0.5 * (M1 - M2)
               mu_range <- c(mu, mu)
             } else{
               mu <- 0.5 * (M1 - M2)
               mu_range <- c(mu, Inf)
             }
           },
           "elemax" = {
             
             hran_diag <- (max(diag(A)) - min(diag(A))) * 0.5
             max_off <- max(abs(A)[row(A)!=col(A)])
             
             if(hran_diag > max_off){
               mu <- 0.5 * (max(diag(A)) + min(diag(A)))
             } else{
               mu <- min(diag(A)) + max_off
             }
             mu_range <- c(mu, mu)
           })
  } else{
    mu_range <- "mu is given by the user"
  }
  alpha <- (mu - lwr_bd) / (mu - res_EV$smallest)
  if(alpha > 1 | alpha < 0){
    alpha <- NA
    message(sprintf("lwr_bd is too large. Should be <= mu (=%.3f)", mu))
  }
  const <- 1
  LPD <- alpha * A + 
    const * mu * (1 - alpha) * diag(ncol(A))
  
  return(list(alpha = alpha,
              mu = mu,
              mu_range = mu_range,
              LPD = LPD))
}

## Depends on: "RSpectra"
## https://cran.r-project.org/web/packages/RSpectra/vignettes/introduction.html
## Find extreme eigenvalues (smallest, largest)
##
## input : (1) A: symmetric matrix
getExtremeEV <- function(A){
  M <- eigs_sym(A, 1, opts = list(retvec = FALSE))$values # largest eigenvalue in absolute value
  m <- eigs_sym(A - M * diag(1, nrow(A)), 1, opts = list(retvec = FALSE))$values
  
  return(list(largest = max(M, M + m),
              smallest = min(M, M + m))
  )
}
## examples
# getExtremeEV(diag(c(10,1,1)))
# getExtremeEV(diag(c(10,1,-2)))
# getExtremeEV(diag(c(3,1,-10)))
# getExtremeEV(diag(c(-2, -5,-10)))
