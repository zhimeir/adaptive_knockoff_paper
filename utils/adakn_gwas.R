adakn_gwas <-  function(W,z,alpha =0.1,offset=1,perc = 0.5,weights = c(1,1),df_cand = 1:5){
  p = length(W)
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  tau.sel = c()
  acc = c()
  
  W_abs = abs(W)
  W_sign = as.numeric(W>0)
  revealed_sign = rep(1,p)
  tau = rep(0,p)
  all_id = 1:p
  
  #initialization
  ref_val = cbind(-z/sum(z),abs(W)/sum(abs(W)))%*%weights
  revealed_id = which(ref_val<=quantile(ref_val,perc))
  print(length(revealed_id))
  #revealed_id = which(W_abs<=s0)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  tau[revealed_id] = W_abs[revealed_id]+1
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  aic = Inf
  df = df_cand[1]
  for(df_id in df_cand){
    mdl  =gam(revealed_sign~ ns(z,df_id) + abs(W),family = "binomial")
    if(df_id - logLik(mdl) <aic){
      aic = df_id - logLik(mdl)
      df = df_id
    }
  }
  print(df)
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl  =gam(revealed_sign~ ns(z,df) + abs(W),family = "binomial")
      fitted.pval = mdl$fitted.values
      plot(fitted.pval)
      fitted.pval = fitted.pval[unrevealed_id]
      predicted.sign = fitted.pval>0.5
      acc  = c(acc, sum(predicted.sign == W_sign[unrevealed_id])/length(unrevealed_id))
      
      ind.min = which(fitted.pval == min(fitted.pval))
      if(length(ind.min)==1){
        ind.reveal = ind.min
      }else{
        ind.reveal = ind.min[which.min(W_abs[ind.min])]
      }
      
      selected.pval = c(selected.pval,min(fitted.pval))
      ind.reveal = unrevealed_id[ind.reveal]
      revealed_id = c(revealed_id,ind.reveal)
      rej.path = c(rej.path,ind.reveal)
      unrevealed_id = all_id[-revealed_id]
      revealed_sign[ind.reveal] = W_sign[ind.reveal]
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      print(fdphat)
      if(fdphat<=fdr){break}
    }
    
    rej = which(W>=tau)  
    rejs[[talpha]] = rej
    nrejs[[talpha]] = length(rej)
    tau.sel = c(tau.sel,length(revealed_id))
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,sel.prob = selected.pval,tau = tau.sel,acc = acc)
  return(result) 
}
