#Soft thresholding estimator
#input--------------------------------------------------------------------------
#data   : given data

#output-------------------------------------------------------------------------
#mat    : matrix with threshold
#idxmat : index matrix with non zero entry


#install.packages("caret")
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
