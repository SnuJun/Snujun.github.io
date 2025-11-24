# evaluation : (A2) Frobenius norm assumption ----------------------------------

eval_ftn_1 <- function(cov.true, cov.est) {
  # (1) Eigenvalues of true covariance matrix
  eig_tr <- eigen(cov.true)$values
  eig_tr <- sort(eig_tr, decreasing = TRUE)
  eig_tr_min <- min(eig_tr)
  eig_tr_m <- mean(eig_tr)
  
  # (2) Eigenvalues of estimated covariance matrix
  eig_est <- eigen(cov.est)$values
  eig_est <- sort(eig_est, decreasing = TRUE)
  eig_est_min <- min(eig_est)
  eig_est_m <- mean(eig_est)
  
  # (3) Comparison
  eig_diff_min <- abs(eig_tr_min - eig_est_min)
  eig_diff_sd <- (nrow(cov.true) - 1) / nrow(cov.true) * sd(eig_tr - eig_est)
  ratio <- eig_diff_min / eig_diff_sd
  
  # Return as named list
  return(list(
    eig_tr_min = eig_tr_min,
    eig_tr_mean = eig_tr_m,
    eig_est_min = eig_est_min,
    eig_est_mean = eig_est_m,
    eig_diff_min = eig_diff_min,
    eig_diff_sd = eig_diff_sd,
    ratio = ratio
  ))
}

# evaluation : Optimal rate over elementwise max norm --------------------------

eval_ftn_2 <- function(cov.true, cov.est, cov.LPD) {
  nom <- max(abs(cov.LPD - cov.true))     # 분자: LPD 오차
  denom <- max(abs(cov.est - cov.true))  # 분모: 기존 오차
  return(nom / denom)
}
