#library(parallel)
library(dplyr)
# basic settings
source("./src_R/datagenerate.R")        # data generating
source("./src_R/threshold.R")           # thresholding estimator
source("./src_R/get_LPD.R")             # LPD function

# evaluation1 ------------------------------------------------------------------

note1 <- function(est){
  eps=1e-4
  eigen_min <- min(eigen(est)$values)
  is_pd <- (eigen_min < eps)
  return(c(eigen_min, is_pd))
}

tot_ev1 <- function(iter, m, n, p) {
  eps=1e-4
  cov.true <- do.call(paste("cov.mtx", m, sep = "."), list(p))
  
  eig.cov.true <- eigen(cov.true)
  sqrt.cov.true <- eig.cov.true$vectors %*% diag(sqrt(eig.cov.true$values)) %*% t(eig.cov.true$vectors)
  
  # sequential processing
  ev_list <- lapply(1:iter, function(i) {
    data.gen <- matrix(rnorm(n*p), nrow=n, ncol=p) %*% sqrt.cov.true
    data <- apply(data.gen, 2, scale, center=TRUE, scale=FALSE)
    cov.soft <- hard_thresholding_CV(data)$mat
    note1(cov.soft)
  })
  
  ev <- do.call(rbind, ev_list)  
  col_means <- mean(ev[,1])
  col_sds <- sd(ev[,1])
  col_sum_div_100 <- sum(ev[,2]) / iter
  
  res <- data.frame(avg_eig_min = col_means, sd_eig_min = col_sds, prop_of_non_PD = col_sum_div_100)
  print(paste(m, '_n=', n, '_p=', p, '_hard thresholding', sep = ""))
  return(res)
}





# evaluation2 ------------------------------------------------------------------
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

tot_ev2 <- function(iter,estim,m,n,p,eps=1e-4){
  set.seed(1234)  #crucial for setting as same model & sample
  cov.true<- do.call(paste("cov.mtx", m, sep = "."), list(p))
  
  eig.cov.true = eigen(cov.true) ;
  eigvec.cov.true = eig.cov.true$vectors ;
  eigval.cov.true = eig.cov.true$values ;
  sqrt.cov.true = eigvec.cov.true %*% diag(sqrt(eigval.cov.true)) %*% t(eigvec.cov.true) ;
  maxeigval.cov.true = max(eigval.cov.true) ;
  mineigval.cov.true = min(eigval.cov.true) ;
  
  if(estim=="Hard"){
    print("implementing with Hard threshold")
  }else if(estim =="Soft"){
    print("implementing with Soft threshold")
  }else if(estim =="Xue"){
    print("implementing with Xue method")
  }else if(estim=="Roth"){
    print("implementing with Rothman method")
  }else if(estim=="LPD_2"){
    print("implementing with LPD_2")
  }else if(estim=="LPD_F"){
    print("implementing with LPD_F")
  }else if(estim == "LPD_1"){
    print("implementing with LPD_1")
  }else if(estim == "LPD_inf"){
    print("implementing with LPD_inf")
  }else if(estim == "COMET"){
    print("implementing with COMET")
  }else{warning("Undefined estimator")}
  
  ev=NULL
  
  for(i in 1:iter){
    cat(m,"_",n,"_",p,"_",estim," with iter : ",i,"\n",sep="")
    data.gen = matrix(rnorm(n*p), nrow=n, ncol=p) ;
    data.gen = data.gen %*% sqrt.cov.true ;
    ## centering data
    data = apply(data.gen, 2, scale, center=TRUE, scale=FALSE) ;
    if(estim=="Hard"){
      est <- hard_thresholding_CV(data)$mat
    }else if(estim =="Soft"){
      est <- soft_thresholding_CV(data)$mat
    }else if(estim =="Xue"){
      est <- ftn.cov.Xue(data,thres.soft,seq(0,1,length=100),eps)$est
    }else if(estim=="Roth"){
      est <- ftn.cov.Roth(data,seq(0,1,length=100))$est
    }else if(estim=="LPD_2"){
      est <- hard_thresholding_CV(data)$mat
      est <- getLPD(est,eps,"spectral")$LPD
    }else if(estim=="LPD_F"){
      est <- hard_thresholding_CV(data)$mat
      est <- getLPD(est,eps,"Frobenius")$LPD
    }else if(estim == "LPD_1"){
      est <- hard_thresholding_CV(data)$mat
      est <- getLPD(est,eps,"elemax")$LPD
    }else if(estim == "LPD_inf"){
      est <- hard_thresholding_CV(data)$mat
      est <- getLPD(est,eps,"Linf")$LPD
    }else if(estim == "COMET"){
      amat <- ftn.cov.adap(data,thresseq = seq(0,1,length=100))$cutmat
      est <- ICF(amat,data)$mat
    }else{warning("Undefined estimator")}

    
    ev=rbind(ev,note2(est,cov.true))
  }
  col_means=colMeans(ev[,1:7])
  col_sds=apply(ev[,1:7],2,sd)
  col_sum_div_100=sum(ev[,8])/iter
  
  res=c(col_means,col_sds,col_sum_div_100)
  res=as.data.frame(t(res))
  colnames(res) <- c('avg_F','avg_S','avg_max','avg_inf','avg_FP','avg_TP','avg_Neigen','sd_F','sd_S','sd_max','sd_inf','sd_FP','sd_TP','sd_Neigen','Pd')
  print(paste(m,'_n=',n,'_p=',p,sep = ""))
  return(res)
}








