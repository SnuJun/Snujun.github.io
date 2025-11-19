#total evaluation2
#This page we'll compare LPD and Comet in several distances

#Base---------------------------------------------------------------------------
library(MASS)     #ginv
library(caret)    #CV
library(RSpectra) #LPD
library(Matrix)   #cov mtx

#Covariance matrix--------------------------------------------------------------

# Model 1: sigma_ij = (1 - |i - j|/10)+
cov.mtx.m1 <- function(p) {
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- max(0, 1 - abs(i - j) / 10)
    }
  }
  return(Sigma)
}
# Model 2: Block matrix
cov.mtx.m2 <- function(p) {
  K <- p / 10  
  if (K != floor(K)) {
    stop("p must be splitted with 10")
  }
  subset_size <- 10
  I_k_max <- seq(10, p, by = 10)  
  
  Sigma <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      indicator1 <- sum(sapply(I_k_max, function(max_idx) {
        (i <= max_idx && j <= max_idx) && ((i > max_idx - subset_size) && (j > max_idx - subset_size))
      }))
      indicator2 <- sum(sapply(I_k_max[-length(I_k_max)], function(max_idx) {
        (i == max_idx + 1 && j == max_idx + 1)
      }))
      
      Sigma[i, j] <- 0.6 * (i == j) + 0.4 * indicator1 + 0.4 * indicator2
    }
  }
  
  return(Sigma)
}
# Model 3: AR(rho=0.7)
cov.mtx.m3 <- function(p){
  eps=1e-3
  mat <- matrix(0, nrow = p, ncol = p)
  diag(mat) <- 1
  for(i in 1:(p-1)){
    mat[i,(i+1):p]=0.7*rbinom(p-i,1,0.01)
    mat[(i+1):p,i]=mat[i,(i+1):p]
  }
  while (min(eigen(mat)$value)<eps) {
    mat=mat+0.01*diag(p)
  }
  return(mat)
}
# Model 4: Sparse model(B(0.3),rho=0.7)
cov.mtx.m4 <- function(p){
  eps=1e-3
  mat <- matrix(0, nrow = p, ncol = p)
  diag(mat) <- 1
  #diag(mat) <- rchisq(p,1)
  for(i in 1:(p-1)){
    mat[i,(i+1):p]=0.7*rbinom(p-i,1,0.3)
    mat[(i+1):p,i]=mat[i,(i+1):p]
  }
  while (min(eigen(mat)$value)<eps) {
    mat=mat+0.01*diag(p)
  }
  return(mat)
}

#functions(Threshold functions)-------------------------------------------------
#(0)

thres.hard = function(obj, cut) {
  ## Note that the threshold can varies upon each element
  res = obj ;
  res[abs(res) < cut] = 0 ;
  return(res) ;
}

thres.soft = function(obj, cut) {
  ## Note that the threshold can varies upon each element
  res = abs(obj) - cut ;
  res[res < 0] = 0 ;
  res = res*sign(obj) ;
  return(res) ;
}

ftn.sparsity = function(symmat) {
  return( sum(abs(symmat) > 10^-6) / length(symmat) ) ;
}

#(1)Hard threshold
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

#(2)Soft threshold
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

#(3)Adaptive hard threshold
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
  attr(resmat, 'cutmat') <- cutmat ;
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

#functions(LPD)-----------------------------------------------------------------
getLPD <- function(A, lwr_bd, type_norm, mu = NULL){
  
  res_EV <- getExtremeEV(A)
  
  if(res_EV$smallest > lwr_bd){
    #message("Smallest eigenvalue is larger than lwr_bd. Returning A unchanged.")
    return(list(alpha = 1,
                mu = NA,
                mu_range = NA,
                LPD = A))
  }
  
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
               #warning(sprintf("The loss function is decreasing in mu.\nA very large mu %.1e is returned.", mu))
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

getExtremeEV <- function(A){
  M <- eigs_sym(A, 1, opts = list(retvec = FALSE))$values # largest eigenvalue in absolute value
  m <- eigs_sym(A - M * diag(1, nrow(A)), 1, opts = list(retvec = FALSE))$values
  
  return(list(largest = max(M, M + m),
              smallest = min(M, M + m))
  )
}

#functions(Comet)---------------------------------------------------------------
ICF <- function(amat, dat, cc=1e-4, tol=1e-04, MAXiter=1e+4,n_type='F'){
  if(sum(amat)==0){
    mat=diag(apply(dat, 2, var))
    return(list('mat'=mat,'iter'=0,'m_df'=0))
  }
  n = nrow(dat); p = ncol(dat)
  
  
  #set initial estimator(zero entry + pd)
  #candidate1 : diagonal mtx
  mat = diag(diag(cov(dat)))
  
  # #candidate2 : sample cov mtx with pdness
  # S = ((n-1)/n)*stats::cov(dat)
  # 
  # mat = S+(cc+abs(min(eigen(mat)$values)))*diag(p)
  
  
  #setting tol and iteration
  matdiff = 1
  iter=1
  #&iter<=100
  while (matdiff > tol&iter <= MAXiter){
    #cat("total iteration:", iter, "\n")
    mat_old = mat
    for (i in 1:p){
      #B = mat[-i,i]*0
      #z = dat[,-i] %*% solve(mat[-i,-i])
      #z = demean(dat[,-i]) %*% solve(mat[-i,-i])
      
      if (sum(amat[i,-i])!=0){
        
        sp_i=drop(which(amat[i,-i]!=0))
        P = as.matrix(ginv(mat[-i,-i])[,sp_i])
        z=dat[,-i]%*%P
        j=crossprod(z,z)/n
        
        #to make pd
        coc=max(abs(eigen(cov(dat))$values[p-1]),cc)
        test=min(eigen(j+coc*crossprod(P,P))$values)
        while(test<cc){
          coc=coc*2
          test=min(eigen(j+coc*crossprod(P,P))$values)
        }
        #cat("var:", i,"eps:",coc, "\n")
        
        off_d=solve(j+coc*crossprod(P,P))%*%t(P)%*%t(dat[,-i])%*%dat[,i]/n
        
        off_dd=as.matrix(ginv(mat[-i,-i])[sp_i,sp_i])
        d=dat[,i]-dat[,-i]%*%P%*%off_d;d=crossprod(d,d)/n
        d=d+coc*(1+crossprod(P%*%off_d)) + t(off_d)%*%off_dd%*%off_d
        
        B = mat[-i,i]*0
        B[as.numeric(which(amat[-i,i]!=0))] = off_d
        B[as.numeric(which(amat[-i,i]==0))] = 0
        mat[i,i] = drop(d)
        mat[i,-i] = B;mat[-i,i] = t(mat[i,-i])
      }
      
    }
    e=min(eigen(mat)$values)
    #cat(iter,"th iteration smallest eigenvalue:",e, "\n")
    
    matdiff =norm(mat-mat_old,type=n_type)
    
    #sum(abs(mat-mat_old))
    #print(matdiff)
    #print(iter)
    iter=iter+1
  }
  return(list('mat'=mat,'iter'=iter,'m_df'=matdiff))
}


#Algorithm (evaluation) --------------------------------------------------------
note2 <- function(est,true){
  tau=1e-5
  res=NULL
  for(t in c("F","2","m","I")){
    # empirical error
    res=c(res,norm(est-true,type=t))
  }
  
  #False positive rate
  res=c(res,sum(abs(est)>=tau & abs(true)<tau) / sum(abs(true)<tau))
  #True positive rate
  res=c(res,sum(abs(est)<tau & abs(true)>=tau) / sum(abs(true)>=tau))
  
  ev=eigen(est)$values
  negeig=sum(abs(ev[ev<tau]))
  #Neg eigenvalues
  res=c(res,negeig)
  #Positive definiteness
  res=c(res,negeig==0)
  
  return(res)
  
}

tot_ev2 <- function(m, n, p, thres = "Hard", iter = 100, n_type = "F", eps = 1e-4) {
  set.seed(1234)
  cov.true <- do.call(paste("cov.mtx", m, sep = "."), list(p))
  eig <- eigen(cov.true)
  sqrt.cov.true <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  
  est_set <- c("init", "LPD_F", "LPD_S", "LPD_inf", "LPD_max", "Comet")
  results <- setNames(vector("list", length(est_set)), est_set)
  com_res <- NULL
  
  for (i in 1:iter) {
    cat(m, "_", n, "_", p, "_with iter: ", i, "\n", sep = "")
    data.gen <- matrix(rnorm(n * p), nrow = n) %*% sqrt.cov.true
    data <- scale(data.gen, center = TRUE, scale = FALSE)
    
    # Threshold 방식 선택
    if (thres == "Hard") {
      ht <- hard_thresholding_CV(data, n_type)
    } else if (thres == "Soft") {
      ht <- soft_thresholding_CV(data, n_type)
    } else {
      ht <- list(
        mat = ftn.cov.adap(data, thresseq = seq(0, 1, length = 100))$estimate,
        idxmat = NULL
      )
      ht$idxmat <- (ht$mat != 0) * 1 - diag(p)
    }
    
    init <- ht$mat
    init_amat <- ht$idxmat
    
    comet_out <- ICF(init_amat,data,n_type=n_type)
    
    estimates <- list(
      init = init,
      LPD_F = getLPD(init, eps, "Frobenius")$LPD,
      LPD_S = getLPD(init, eps, "spectral")$LPD,
      LPD_inf = getLPD(init, eps, "Linf")$LPD,
      LPD_max = getLPD(init, eps, "elemax")$LPD,
      Comet = comet_out$mat
    )
    
    com_res <- rbind(com_res, c(comet_out$iter, comet_out$m_df))
    
    for (estim in est_set) {
      results[[estim]] <- rbind(results[[estim]], note2(estimates[[estim]], cov.true))
    }
  }
  
  
  tot_res <- NULL
  for (estim in est_set) {
    ev <- results[[estim]]
    col_means <- colMeans(ev[, 1:7])
    col_sds <- apply(ev[, 1:7], 2, sd)
    pd_rate <- sum(ev[, 8]) / iter
    res <- data.frame(t(c(col_means, col_sds, pd_rate)))
    colnames(res) <- c('avg_F', 'avg_S', 'avg_max', 'avg_inf', 'avg_FP', 'avg_TP', 'avg_Neigen',
                       'sd_F', 'sd_S', 'sd_max', 'sd_inf', 'sd_FP', 'sd_TP', 'sd_Neigen', 'Pd')
    rownames(res) <- estim
    tot_res <- rbind(tot_res, res)
  }
  
  cat("Total result w.r.t ", n_type, " norm\n", sep = "")
  return(list('tot_res'=tot_res,"com_res"=com_res))
}


#main---------------------------------------------------------------------------
thres="Hard"
set.seed(99)
n = 100
p_set = c(30,60,120)
dist_set=c("F","2","m","I")
iter = 100
m_set = c('m1','m2','m3','m4')

result=NULL
for(m in m_set){
  for(p in p_set){
    for(dist in dist_set){
      result=rbind(result,tot_ev2(m,n,p,thres=thres,iter=iter,n_type=dist))
    }
  }
}
print(result)
write.csv(result, "totresult.csv", row.names = FALSE)