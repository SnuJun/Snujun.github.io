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
  
  ## If A already satisfies the lower eigenvalue bound, return A itself
  if (res_EV$smallest >= lwr_bd) {
    return(list(
      alpha = 1,
      mu = NA,
      mu_range = NA,
      LPD = A,
      status = "already_PD"
    ))
  }
  
  if(is.null(mu)){
    
    switch(as.character(type_norm),
           
           "spectral" = {
             mu <- max(
               0.5 * (res_EV$largest + res_EV$smallest),
               lwr_bd
             )
             mu_range <- c(mu, Inf)
           },
           
           "Frobenius" = {
             res_eig <- eigen(A, only.values = TRUE)
             lam1 <- min(res_eig$values)
             
             mu <- sum((res_eig$values - lam1)^2) /
               sum((res_eig$values - lam1))
             
             mu <- max(mu, lwr_bd)
             mu_range <- c(mu, mu)
           },
           
           "Linf" = {
             M1 <- max(rowSums(abs(A)))
             
             A_tmp <- A - 2 * diag(diag(A))
             A_tmp[row(A_tmp) != col(A_tmp)] <-
               abs(A_tmp[row(A_tmp) != col(A_tmp)])
             
             M2 <- max(rowSums(A_tmp))
             
             if(res_EV$smallest + M2 > 1e-10){
               
               mu <- Inf
               mu_range <- c(Inf, Inf)
               
               LPD <- A + (lwr_bd - res_EV$smallest) * diag(ncol(A))
               
               return(list(
                 alpha = 1,
                 mu = mu,
                 mu_range = mu_range,
                 LPD = LPD,
                 status = "mu_infinity_diagonal_shift"
               ))
               
             } else if(res_EV$smallest + M2 < -1e-10){
               
               mu_raw <- 0.5 * (M1 - M2)
               mu <- max(lwr_bd, mu_raw)
               mu_range <- c(mu, mu)
               
             } else{
               
               mu_raw <- 0.5 * (M1 - M2)
               mu <- max(lwr_bd, mu_raw)
               mu_range <- c(mu, Inf)
             }
           },
           
           "elemax" = {
             
             hran_diag <- (max(diag(A)) - min(diag(A))) * 0.5
             max_off <- max(abs(A)[row(A) != col(A)])
             
             if(hran_diag > max_off){
               mu <- 0.5 * (max(diag(A)) + min(diag(A)))
             } else{
               mu <- min(diag(A)) + max_off
             }
             
             mu <- max(mu, lwr_bd)
             mu_range <- c(mu, mu)
           })
    
  } else{
    mu_range <- "mu is given by the user"
  }
  
  if(!is.finite(mu)){
    
    LPD <- A + (lwr_bd - res_EV$smallest) * diag(ncol(A))
    
    return(list(
      alpha = 1,
      mu = mu,
      mu_range = mu_range,
      LPD = LPD,
      status = "mu_infinity_diagonal_shift"
    ))
  }
  
  alpha <- (mu - lwr_bd) / (mu - res_EV$smallest)
  
  if(is.na(alpha) || alpha > 1 || alpha < 0){
    stop(sprintf(
      "Invalid alpha: alpha = %.6f, mu = %.6f, lwr_bd = %.6f, lambda_min = %.6f",
      alpha, mu, lwr_bd, res_EV$smallest
    ))
  }
  
  LPD <- alpha * A + mu * (1 - alpha) * diag(ncol(A))
  
  return(list(
    alpha = alpha,
    mu = mu,
    mu_range = mu_range,
    LPD = LPD,
    status = "success"
  ))
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
