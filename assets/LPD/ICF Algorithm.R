#ICF algorithm with p>=n
library(MASS)   # ginv

#input--------------------------------------------------------------------------
#cc     : constants that make initial estimator(sample covariance) pd
#amat   : idx matrix which is nonzero(defined by Hard threshold)
#dat    : given data
#tol    : constants for stoping iteration

#output-------------------------------------------------------------------------
#mat    : matrix with ICF algorithm



#deriving non zero mat with adaptive hard threshold ----------------------------
ftn.sparsity = function(symmat) {
  return( sum(abs(symmat) > 10^-6) / length(symmat) ) ;
}


thres.hard = function(obj, cut) {
  ## Note that the threshold can varies upon each element
  res = obj ;
  res[abs(res) < cut] = 0 ;
  return(res) ;
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
  attr(resmat, 'cutmat') <- cutmat ;
  return(resmat) ;
}

ftn.cov.adap = function(datamat, thres.ftn= thres.hard,
                        thresseq = NULL, thresnum = 100, diag.thres=FALSE) {
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


# ICF algrotithm with general setting ------------------------------------------
ICF <- function(amat, dat, cc=1e-4, tol=1e-02, MAXiter=1e+3){
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
    matdiff =norm(mat-mat_old,'F')
    
    #sum(abs(mat-mat_old))
    #print(matdiff)
    #print(iter)
    iter=iter+1
  }
  return(list('mat'=mat,'iter'=iter,'m_df'=matdiff))
}


# demean = function(dat){
#   meanmat = matrix(rep(colMeans(dat), nrow(dat)), ncol = ncol(dat), byrow = TRUE)
#   datdm = dat - meanmat
#   return(datdm)
# }



