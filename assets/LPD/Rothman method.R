#Rothman method

library(PDSCE)

ftn.cov.Roth = function(datamat, thresseq, tol = 10^-7, MAXITER = 10000) {
  time1 = Sys.time() ;
  roth = pdsoft.cv(x = datamat, lam.vec = thresseq, standard = FALSE,init="soft", nsplits = 5, n.tr = 4, tolout = tol, maxitout=MAXITER,quiet=TRUE) ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(est = roth$sigma, iter = NA, tol = tol,
              thres = roth$best.lam, time = time2)) ;
}


