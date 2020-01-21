### set up all methods
# adapted from https://github.com/lihualei71/adaptPaper/blob/master/R/other_methods.R

BH_method = function(pvals,alpha){
  khat=max(c(0,which(sort(pvals)<=alpha*(1:length(pvals))/length(pvals))))
  which(pvals<=alpha*khat/length(pvals))
}

Storey_method = function(pvals,alpha,thr){
  est_proportion_nulls=min(1,mean(pvals>thr)/(1-thr))
  pvals[pvals>thr] = Inf
  khat=max(c(0,which(sort(pvals)<=alpha/est_proportion_nulls*(1:length(pvals))/length(pvals))))
  which(pvals<=alpha/est_proportion_nulls*khat/length(pvals))
}


SABHA_method = function(pvals, qhat, alpha, tau){
  # Use the original, or estimated q as input
  pvals[pvals>tau] = Inf
  khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
  which(qhat*pvals<=alpha*khat/length(pvals))
}

Adaptive_SeqStep_method = function(pvals, alpha, thr1, thr2){ # Lei & Fithian 2016's method
  # thr1 & thr2 correspond to s & lambda (with s<=lambda) in their paper
  fdphat = thr1 / (1-thr2) * (1 + cumsum(pvals>thr2)) / pmax(1, cumsum(pvals<=thr1))
  if(any(fdphat<=alpha)){
    khat = max(which(fdphat<=alpha))
    return(which(pvals[1:khat]<=thr1))
  }else{
    return(NULL)
  }

}

