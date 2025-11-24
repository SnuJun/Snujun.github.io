# basic settings
source("./simulation/src_R/datagenerate.R")       # data generating
source("./simulation/src_R/threshold.R")          # thresholding estimator
source('/./simulation/evaluation.R')              # evaluation
source("./src_R/get_LPD.R")             # LPD function

# simulation 1 -----------------------------------------------------------------

set.seed(1234)
n = 100
p_set = c(30,60,120)
iter = 100
m_set = c( 'm1','m2','m3','m4') 
result=NULL
for (m in m_set) {
  datum <- do.call(rbind, lapply(p_set, function(p) tot_ev1(iter, m, n, p)))
  result=rbind(result,datum)
}


# simulation 2 -----------------------------------------------------------------

set.seed(99)
n = 100
p_set = c(30,60,120)
iter = 100
m_set = c('m1','m2','m3', 'm4')
estim_set=c("Hard","Soft","LPD_F","LPD_2","LPD_1","LPD_inf")
#"Roth","Xue","COMET"
result=NULL
for(m in m_set){
  for(p in p_set){
    for(estim in estim_set){
      result=rbind(result,tot_ev2(iter,estim,m,n,p))
    }
  }
}

