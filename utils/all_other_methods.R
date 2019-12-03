### set up all methods


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


#num_alpha=length(alphalist); max_alpha = max(alphalist)

# gather results
#NumRej = array(0,c(7,num_alpha,3))
# methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 SABHA
#for(j in 1:3){
#  if(j==1){pvals=output$pvals[output$ord]}else{if(j==2){pvals=output$pvals[output$ord_small]}else{pvals=output$pvals[output$ord_random]}}
#  qhat = Solve_q_step(pvals,tau,eps)
#  for(i in 1:num_alpha){
#    NumRej[1,i,j] = SeqStep(pvals,alpha=alphalist[i],C=2)
#    NumRej[2,i,j] = HingeExp(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i],C=2)
#    NumRej[3,i,j] = ForwardStop(pvals*(1-1/(1+choose(10,5))),alpha=alphalist[i])
#    NumRej[4,i,j] = length(Adaptive_SeqStep_method(pvals, alphalist[i], thr1, thr2))
#    NumRej[5,i,j] = length(BH_method(pvals, alphalist[i]))
#    NumRej[6,i,j] = length(Storey_method(pvals, alphalist[i], thr))
#    NumRej[7,i,j] = length(SABHA_method(pvals, qhat, alphalist[i], tau))
#  }
#}

#gene_drug_plot_results(NumRej,'gene_drug_results.pdf')
