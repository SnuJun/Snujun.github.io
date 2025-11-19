#Xue method
# Step1 : deriving soft threshold matrix 
ftn.posi.mat = function(mat, mineigval) {
  ### Projection onto the PSD cone
  #  if(!isSymmetric(mat)) {
  #    if(isSymmetric(mat, tol=10^-3)) { mat = (mat + t(mat)) / 2 ; } else {
  #    stop("Input is not symmetric")  ; }
  #  }
  
  eig = eigen(mat) ;
  eigval = eig$values ;
  if (sum(eigval < mineigval) == 0) return(list(res=mat, operated=c(FALSE))) ;
  eigval[eigval < mineigval] = mineigval ;
  eigvec = eig$vectors ;
  
  res = eigvec %*% diag(eigval) %*% t(eigvec) ;
  res = (res + t(res)) / 2 ;  # to stabilize symmetricity
  return(list(res=res, operated=c(TRUE))) ; 
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


subftn.unithres = function(datamat=NULL, sampcov, thres.ftn=thres.soft, cut, diag.thres=FALSE) {
  
  resmat = thres.ftn(obj = sampcov, cut = cut) ;
  if (!diag.thres) diag(resmat) = diag(sampcov) ;
  return(resmat) ;
}


ftn.cov.uni = function(datamat, thres.ftn=thres.soft, thresseq=NULL,thresnum=100, diag.thres=FALSE) {
  #thresseq have priority
  
  ## (cutpoint estimation)
  ## split the samples
  n = nrow(datamat) ; p = ncol(datamat) ;
  ## centerizing
  data = apply(datamat, 2, scale, center=T, scale=F) ;
  cov.samp = t(data) %*% data / (n-1) ;
  cov.samp.off = cov.samp ; diag(cov.samp.off) = 0 ;
  max.thres = max(abs(cov.samp.off)) ;
  thresseq.temp = 1.05^seq(from=0, to=thresnum, by=1) - 1 ;
  if (is.null(thresseq)) thresseq = max.thres * thresseq.temp/max(thresseq.temp) ;
  
  time1 = Sys.time() ;
  n.test = floor(n / 5) ;
  n.fold = 5 ;
  #  J = floor(1 / sqrt(log(p)/n)) ;
  ## make the stack of errors
  #  thresseq = sqrt(log(p)/n) * (0:J) ;
  
  errorseq = rep(0, length(thresseq)) ;
  nzseq = rep(0, length(thresseq)) ;
  ### for each training/test set
  ### randomly choose the validation sets
  set.seed(1) ;
  perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;
  
  for (i in 1:n.fold) {
    test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
    test.ind.perm = perm.ind[test.ind] ;
    data.test = data[test.ind.perm, ] ;
    data.train = data[-test.ind.perm, ] ;
    ### test sample cov matrix and training sample cov matrix
    sampcov.test = var(data.test) ;
    sampcov.train = var(data.train) ;
    
    #### for each candidate : grid on the integer multiplication of log(p)/n
    for ( j in 1:length(thresseq) ) {
      #### training (thresholded) estimator
      #function(datamat, sampcov, thres.ftn, cut, diag.thres=FALSE)
      cov.train = subftn.unithres(datamat=NULL, sampcov=sampcov.train, thres.ftn=thres.ftn,
                                  cut=thresseq[j], diag.thres=diag.thres) ;
      #### frobenius norm loss
      error = sum((cov.train - sampcov.test)^2) ;
      errorseq[j] = errorseq[j] + error ;
      nzseq[j] = nzseq[j] + ftn.sparsity(cov.train) ;
    }
  }
  #time2 = Sys.time() ;
  #print(time2 - time1) ;
  ## (select the optimal cutpoint)
  nzseq = nzseq / n.fold ;
  thres.opt = thresseq[which.min(errorseq)] ;
  cutmat = matrix(thres.opt, p, p) ;
  if (!diag.thres) diag(cutmat) = 0 ;
  
  ## (store the thresholding estimator)
  resmat = subftn.unithres(datamat=NULL, sampcov=cov.samp, thres.ftn=thres.ftn,
                           cut=thres.opt, diag.thres=diag.thres) ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(estimate = resmat, lambda.opt = thres.opt, cutmat = cutmat,
              nz.opt = ftn.sparsity(resmat), lambdaseq=thresseq, errorseq = errorseq,
              nzseq = nzseq, time = time2, iter = NA, tol = NA)) ;
}

# Xue method
subftn.Xue = function(sampcov, initial.cov.thres, cut, mineigval=10^-3, diag.thres=FALSE, tolprimal=10^-5, toldual=10^-4, MAXITER=1000) {
  
  #### MAIN ALGORITHM
  tuningmat.old = 0 ;
  mu = 1 ;
  covmat.old = initial.cov.thres ;
  Soff = ( sum(abs(sampcov)) - sum(abs(diag(sampcov))) ) / 2 ;
  
  errorvec1 = rep(0, MAXITER) #####  
  errorvec2 = rep(0, MAXITER) #####    
  for (iter in 1:MAXITER) {
    Thetamat.new.list = ftn.posi.mat(mat = covmat.old + mu * tuningmat.old, mineigval = mineigval) ;
    # if ((iter == 1) & !(Thetamat.new.list$operated)) {
    #   return(list(est=covmat.old, iter=0, tol=tol)) ; }
    Thetamat.new = Thetamat.new.list$res ;
    
    temp = mu*(sampcov - tuningmat.old) + Thetamat.new ;
    covmat.new = thres.soft( temp, cut*mu ) / (1 + mu) ;
    if (!diag.thres) diag(covmat.new) = diag(temp) / (1 + mu) ;
    
    tuningmat.new = tuningmat.old - (Thetamat.new - covmat.new)/mu ;
    
    #### convergence criteria
    err.normalizer = Soff ;               
    error1 = sum(abs(covmat.new - covmat.old)) / err.normalizer ;
    error2 = sum(abs(covmat.new - Thetamat.new)) / err.normalizer ;
    errorvec1[iter] = error1 ; #######
    errorvec2[iter] = error2 ; #######    
    #if ( error1  < tol) {
    if ( (error1 < tolprimal) & (error2 < toldual) ) {
      #      cat("Funtion subftn.Xue : convergence acheived at step", iter, "\n") ;
      #      windows() ; plot.ts(errorvec1, ylim=c(0,10^-5)) ;
      break ;
    } else {
      covmat.old = covmat.new ;
      tuningmat.old = tuningmat.new ;
    }
    if (iter == MAXITER) { 
      #      cat("Funtion subftn.Xue : reached to maximum iteration\n") ;
      #    windows() ; plot.ts(log10(errorvec1), ylim=c(-6,1)) ;
      #      windows() ; plot.ts(log10(errorvec2), ylim=c(-6,1)) ;
      #      windows() ; plot.ts(errorvec2, ylim=c(0,10^-3)) ;
    }
  }
  
  attr(covmat.new, 'iter') <- iter ;
  attr(covmat.new, 'tol') <- tolprimal ;
  return(covmat.new) ;
}

ftn.cov.Xue = function(datamat, thres.ftn = thres.soft, thresseq, mineigval = 1e-3, tol = 1e-4, MAXITER = 100, message = FALSE) {
  time1 = Sys.time()
  
  n = nrow(datamat)
  p = ncol(datamat)
  n.test = floor(n / 5)
  n.fold = 5
  errorstack = rep(0, length(thresseq))
  
  ## Center data
  data = scale(datamat, center = TRUE, scale = FALSE)
  
  set.seed(1)
  perm.ind = sample.int(n, size = n, replace = FALSE)
  if (message) {
    cat("Permutated indices:\n")
    print(perm.ind)
  }
  
  for (i in 1:n.fold) {
    if (message) cat("Test set: fold", i, "\n")
    
    test.ind = ((i - 1) * n.test + 1):(i * n.test)
    test.ind.perm = perm.ind[test.ind]
    data.test = data[test.ind.perm, ]
    data.train = data[-test.ind.perm, ]
    
    sampcov.test = var(data.test)
    sampcov.train = var(data.train)
    
    for (j in 1:length(thresseq)) {
      if (message) cat("Threshold:", thresseq[j], "\n")
      
      initial.cov.thres = sampcov.train  # could be pre-thresholded, but here use raw
      cov.train = subftn.Xue(sampcov = sampcov.train,initial.cov.thres = initial.cov.thres,cut = thresseq[j],mineigval = mineigval,tolprimal = tol,MAXITER = MAXITER)
      
      error = norm(cov.train - sampcov.test,'F')
      errorstack[j] = errorstack[j] + error
    }
  }
  
  if (message) {
    print(rbind(thresseq, errorstack))
  }
  
  ## Optimal threshold
  thres.opt = thresseq[which.min(errorstack)]
  
  ## Final estimate on full dataset
  cov.samp = var(data)
  cov.Xue = subftn.Xue(sampcov = cov.samp,
                       initial.cov.thres = cov.samp,
                       cut = thres.opt,
                       mineigval = mineigval,
                       tolprimal = tol,
                       MAXITER = MAXITER)
  
  time2 = difftime(Sys.time(), time1, units = "secs")
  return(list(est = cov.Xue,
              thres = thres.opt,
              iter = attr(cov.Xue, "iter"),
              tol = attr(cov.Xue, "tol"),
              time = time2))
}








