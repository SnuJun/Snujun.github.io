#(1) Hard threshold ------------------------------------------------------------
library(caret)
hard_thresholding_CV <- function(data,n_type="F") {
  cov.samp=cov(data)
  cov.diag=diag(diag(cov.samp))
  cov.diag.off=cov.samp-cov.diag
  
  threshold_candidates = seq(0,1,length.out=100)
  
  n <- nrow(data)
  p <- ncol(data)
  
  # 5-Fold Cross Validation
  folds <- createFolds(1:n, k = 5, list = TRUE)
  cv_errors <- numeric(length(threshold_candidates))  
  for (t in seq_along(threshold_candidates)) {  
    thr <- threshold_candidates[t]
    fold_errors <- numeric(length(folds)) 
    for (fold_idx in seq_along(folds)) {  
      test_idx <- folds[[fold_idx]]
      train_idx <- setdiff(1:n, test_idx)
      
      data_train <- data[train_idx, ]
      data_test <- data[test_idx, ]
      
      
      cov.train <- cov(data_train)
      
      
      cov.hard <- cov.train
      cov.hard[abs(cov.hard) < thr] <- 0  
      
      cov.test <- cov(data_test)
      
      error <- norm(cov.hard - cov.test, type = n_type)  # Frobenius norm
      
      fold_errors[fold_idx] <- error  
    }
    
    
    cv_errors[t] <- sum(fold_errors)  
  }
  
  
  best_thr <- threshold_candidates[which.min(cv_errors)]
  #plot(threshold_candidates,cv_errors, main='Case2');abline(v=best_thr)
  #print(paste("Optimal threshold selected by 5-fold CV:", round(best_thr, 4)))
  
  
  cov.final <- cov.diag.off
  idxmat=(abs(cov.final) >= best_thr)*1
  cov.final[abs(cov.final) < best_thr] <- 0;cov.final=cov.final+cov.diag
  return(list('mat'=cov.final,'idxmat'=idxmat))
}

#(2) Soft threshold ------------------------------------------------------------
library(caret)
thres.soft = function(obj, cut) {
  ## Note that the threshold can varies upon each element
  res = abs(obj) - cut ;
  res[res < 0] = 0 ;
  res = res*sign(obj) ;
  return(res) ;
}
soft_thresholding_CV <- function(data,n_type="F") {
  cov.samp=cov(data)
  cov.diag=diag(diag(cov.samp));cov.diag.off=cov.samp-cov.diag
  
  threshold_candidates = seq(0,1,length.out=100)
  
  n <- nrow(data)
  p <- ncol(data)
  
  # 5-Fold Cross Validation
  folds <- createFolds(1:n, k = 5, list = TRUE)
  cv_errors <- numeric(length(threshold_candidates))  
  for (t in seq_along(threshold_candidates)) {  
    thr <- threshold_candidates[t]
    fold_errors <- numeric(length(folds)) 
    for (fold_idx in seq_along(folds)) {  
      test_idx <- folds[[fold_idx]]
      train_idx <- setdiff(1:n, test_idx)
      
      data_train <- data[train_idx, ]
      data_test <- data[test_idx, ]
      
      
      cov.train <- cov(data_train)
      
      cov.soft <- cov.train
      cov.soft.diag <- diag(diag(cov.soft));cov.soft <- cov.soft-cov.soft.diag
      cov.soft <- thres.soft(cov.soft,thr);cov.soft <- cov.soft+cov.soft.diag
      
      cov.test <- cov(data_test)
      
      error <- norm(cov.soft - cov.test, n_type)  
      
      fold_errors[fold_idx] <- error  
    }
    
    
    cv_errors[t] <- sum(fold_errors)  
  }
  
  
  best_thr <- threshold_candidates[which.min(cv_errors)]
  #plot(threshold_candidates,cv_errors, main='Case2');abline(v=best_thr)
  #print(paste("Optimal threshold selected by 5-fold CV:", round(best_thr, 4)))
  
  
  cov.final <- cov.diag.off
  idxmat=(abs(cov.final) >= best_thr)*1
  cov.final[abs(cov.final) < best_thr] <- 0;cov.final=cov.final+cov.diag
  return(list('mat'=cov.final,'idxmat'=idxmat))
}

#(3) Adative Hard threshold ----------------------------------------------------

thres.hard <- function(obj, cut) {
  res <- obj
  res[abs(res) < cut] <- 0
  return(res)
}
ftn.sparsity <- function(symmat) {
  return(sum(abs(symmat) > 1e-6) / length(symmat))
}
subftn.adapthres = function(datamat, sampcov, thres.ftn=thres.hard, cut, diag.thres=FALSE) {
  n = nrow(datamat) ; p = ncol(datamat) ;
  datamat.centered = datamat - rep(colMeans(datamat), rep(n, p)) ;
  datamat.centered = datamat.centered^2 ;
  ## sqrt(theta_ij)
  stdevmat.sampcov = sqrt( t(datamat.centered) %*% datamat.centered / n  - 
                             (n - 2) / n * sampcov^2 ) ;
  #  print(stdevmat.sampcov)
  cutmat = cut * stdevmat.sampcov ;
  if(!diag.thres) diag(cutmat) = 0 ;
  
  #  print(cutmat) ;
  resmat = thres.ftn(obj = sampcov, cut = cutmat) ;
  return(resmat) ;
}
ftn.cov.adap = function(datamat, thres.ftn= thres.hard,thresseq = NULL, thresnum = 100, diag.thres=FALSE) {
  thresseq.temp = 1.05^seq(from=0, to=thresnum, by=1) - 1 ;
  if (is.null(thresseq)) thresseq = 4 * thresseq.temp/max(thresseq.temp) ;
  
  time1 = Sys.time() ;
  
  ## (cutpoint estimation)
  ## split the samples
  n = nrow(datamat) ; p = ncol(datamat) ;
  n.test = floor(n / 5) ;
  n.fold = 5 ;
  
  ## HERE IS DIFFERENT FROM OTHERS
  ## make the stack of errors
  
  errorseq = rep(0, length(thresseq)) ;
  nzseq = rep(0, length(thresseq)) ;
  ### for each training/test set
  ### randomly choose the validation sets
  set.seed(1) ;
  perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;
  
  for (i in 1:n.fold) {
    test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
    test.ind.perm = perm.ind[test.ind] ;
    data.test = datamat[test.ind.perm, ] ;
    data.train = datamat[-test.ind.perm, ] ;
    ### test sample cov matrix and training sample cov matrix
    sampcov.test = var(data.test) ;
    sampcov.train = var(data.train) ;
    
    for ( j in 1:length(thresseq) ) {
      #### training (thresholded) estimator
      cov.train = subftn.adapthres(datamat = data.train, sampcov = sampcov.train,
                                   thres.ftn = thres.ftn, cut = thresseq[j] * sqrt(log(p) / n), diag.thres=diag.thres )  ;
      #### frobenius norm loss
      error = sum((as.numeric(cov.train - sampcov.test))^2) ;
      errorseq[j] = errorseq[j] + error ;
      nzseq[j] = nzseq[j] + ftn.sparsity(cov.train) ;
    }
  }
  #time2 = Sys.time() ;
  #print(time2 - time1) ;
  ## (select the optimal cutpoint)
  nzseq = nzseq / n.fold ;
  thres.opt = thresseq[which.min(errorseq)] ;
  
  
  ## (store the thresholding estimator)
  resmat = subftn.adapthres(datamat = datamat, sampcov = var(datamat),
                            thres.ftn = thres.ftn, cut = thres.opt * sqrt(log(p) / n), diag.thres=diag.thres)  ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(estimate = resmat, lambda.opt = thres.opt, cutmat = attr(resmat, 'cutmat'),
              nz.opt = ftn.sparsity(resmat), lambdaseq=thresseq, errorseq = errorseq,
              nzseq = nzseq, time = time2, iter = NA, tol = NA)) ;                           
}

#(4) Bandable estimator --------------------------------------------------------
ftn.cov.band <- function(X, band.size = 10) {
  n <- nrow(X)
  p <- ncol(X)
  S <- cov(X)
  S_band <- matrix(0, p, p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      if (abs(i - j) <= band.size) {
        S_band[i, j] <- S[i, j]
      }
    }
  }
  
  return(list(estimate = S_band, band.size = band.size))
}

#(5) Blockwise tridiagonal estimator -------------------------------------------
ftn.cov.blck <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  # block size 자동 설정
  block.size <- p / 10
  
  if (block.size != floor(block.size)) {
    stop("p must be divisible by 10.")
  }
  
  K <- 10  # 블록 개수
  S <- cov(X)
  block.indices <- split(1:p, ceiling(seq_along(1:p) / block.size))
  
  S_block <- matrix(0, p, p)
  
  for (i in 1:K) {
    for (j in (i-1):(i+1)) {
      if (j >= 1 && j <= K) {
        rows <- block.indices[[i]]
        cols <- block.indices[[j]]
        S_block[rows, cols] <- S[rows, cols]
      }
    }
  }
  
  return(list(estimate = S_block, block.size = block.size))
}

#(6) POET estimator ------------------------------------------------------------
library(POET)

