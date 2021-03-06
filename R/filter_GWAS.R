filter_GWAS <- function(W,U,alpha = 0.1,offset = 1,mute = TRUE,df =3,R=1,s0=NULL){
  #Check the input format
  if(is.numeric(W)){
    W = as.vector(W)
  }else{
    stop('W is not a numeric vector')
  }
  
  if(is.numeric(U) ==1){
    U = as.matrix(U)
  }else{
    stop('U is not numeric')
  }
  
  #Extract dimensionality
  p = length(W)
  
  #check if z is in the correct form
  if(dim(U)[1]!=p){
    if(dim(U)[2]==p){
      U = t(U)
    }
    else{
      stop('Please check the dimensionality of the side information!')
    }
  }
  
  # dimention of the side information
  pz = dim(U)[2]
  
  # Initializing the output
  all_id = 1:p
  tau.sel = c() #stopping time
  rejs = vector("list",length(alpha)) #rejection sets
  nrejs =  rep(0,length(alpha)) #number of rejections
  ordered_alpha = sort(alpha,decreasing = TRUE) #ordered target levels: from the largest to the smallest
  rej.path = c()
  
  
  # Algorithm initialization
  eps = 1e-10
  W_abs = abs(W)
  W_revealed = W_abs
  tau = W_abs
  revealed_id = which(W_abs<=s0)# Reveal a small amount of signs based on magnitude only
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  W_revealed[revealed_id] = W[revealed_id]
  rej.path = c(rej.path,revealed_id)
  tau[revealed_id] =  W_abs[revealed_id]+1
  }else{
    unrevealed_id = all_id
  }
  
  pi = rep(sum(W>0)/p,p)
  delta0 = sum(W==0)/p*sum(W<=0)/p
  delta1 = sum(W==0)/p*sum(W>0)/p
  t = logis(W_revealed)
  mu_0 =  rep(-log(logis(mean(W_revealed[W_revealed<0]))),p)
  mu_1 = rep(-log(logis(mean(W_revealed[W_revealed>0]))),p)
  mu_1[W==0] = log(2)
  mu_0[W==0] = log(2)
  H = rep(eps,p)
  y0 = -log(t)
  y1 = -log(t)
  
  # compute estimated FDR 
  fdphat = calculate.fdphat(W,tau,offset = offset)
  if(mute == FALSE) print(fdphat)
  
  for (talpha in 1:length(alpha)){
    fdr = ordered_alpha[talpha]
    while(fdphat>fdr && length(unrevealed_id)>=1){
      # EM algorithm
      for (r in 1:R){
        # E step
        # revealed part
        H[revealed_id] = sapply(1:length(revealed_id),function(j) prob_revealed(pi[revealed_id[j]],mu_0[revealed_id[j]],mu_1[revealed_id[j]],t[revealed_id[j]],delta0,delta1))
        y1[revealed_id] = -log(t[revealed_id])
        y0[revealed_id] = -log(t[revealed_id])
        
        # unrevealed part
        H[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) prob_unrevealed(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))
        H = pmax(H,eps) #avoid numerical errors
        H = pmin(H,1-eps) #avoid numerical errors
        y1[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) exp_unrevealed_1(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))/H[unrevealed_id]
        y0[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) exp_unrevealed_0(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))/(1-H[unrevealed_id])
        
        # M step
          #update nu
          if(dim(U)[2] ==1){
            mdl = gam(log(H/(1-H))~ns(U,df))
            pi = logis(mdl$fitted.values)
          }else{
            mdl = gam(log(H/(1-H))~s(U[,1],U[,2]))
            pi = logis(mdl$fitted.values)
          }
          
          #update beta_0
          if(dim(U)[2]==1){mdl =gam(y0[t!=1/2]~ns(U[t!=1/2],df),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
          }else{
            mdl =gam(y0[t!=1/2]~s(U[t!=1/2,1],U[t!=1/2,2]),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
          }
          mu_0[t!=1/2] =mdl$fitted.values
          
          #update beta_1
          if(dim(U)[2]==1){mdl =gam(y1[t!=1/2]~ns(U[t!=1/2],df),weights = (H[t!=1/2]),family = Gamma(link = "log"))
          }else{
            mdl =gam(y1[t!=1/2]~s(U[t!=1/2,1],U[t!=1/2,2]),weights = (H[t!=1/2]),family = Gamma(link = "log"))
          }
          mu_1[t!=1/2] =mdl$fitted.values
        
        #update delta
        delta0 = sum((1-H)*(t==1/2))/(sum((1-H)*(t==1/2))+sum((1-H)*(t!=1/2)))
        delta1 = sum((H)*(t==1/2))/(sum((H)*(t==1/2))+sum((H)*(t!=1/2)))
      }
      
      
      horder = sapply(1:length(unrevealed_id),function(j) order_prob(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))
      ind.min = which(horder == min(horder))
      if(length(ind.min)==1){
        ind.reveal = ind.min
      }else{
        ind.reveal = ind.min[which.min(W_abs[ind.min])]
      }
      
      ind.reveal = unrevealed_id[ind.reveal]
      revealed_id = c(revealed_id,ind.reveal)
      unrevealed_id = all_id[-revealed_id]
      W_revealed[ind.reveal] = W[ind.reveal]
      t = logis(W_revealed)
      rej.path = c(rej.path,ind.reveal)
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      
      # compute estimated FDR
      fdphat = calculate.fdphat(W,tau,offset = offset)
      if(mute == FALSE) print(fdphat)
    }
    
    # determine selection set
      rej = which(W>=tau)
      rejs[[talpha]] = rej
      nrejs[talpha] = length(rej)
      tau.sel = c(tau.sel, length(revealed_id))
    
    
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,tau = tau.sel,W=W)
  return(result)
}


## Utility funcitons
calculate.fdphat = function(W,tau,offset = 1){
  p = length(W)
  fdphat = (offset+sum(W<=-tau))/max(sum(W>=tau),1)
  return(fdphat)
}


logis <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}

dens  <- function(t,mu,delta){
  val = delta*(t==1/2)+(1-delta)*(t!=1/2)*(t^(1/mu-1)/mu)
  return(val)
}


prob_revealed <- function(pi,mu0,mu1,t,delta0,delta1){
  num = pi*dens(t,mu1,delta1)
  denom = pi*dens(t,mu1,delta1) + (1-pi)*dens(t,mu0,delta0)
  if(denom==0){value=0.5}else{value = num/denom}#avoid numerical errors
  return(value)
}

prob_unrevealed <- function(pi,mu0,mu1,t,delta0,delta1){
  num = pi*dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1)
  denom = num +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom==0){value=0.5}else{value = num/denom}#avoid numerical errors
  return(value)
}

exp_unrevealed_0 <- function(pi,mu0,mu1,t,delta0,delta1){
  y1  = -log(t)
  y2 = -log(1-t)
  num = y1 *( 1-pi) *dens(t,mu0,delta0) + y2*(1-pi)*dens(1-t,mu0,delta0)
  denom =pi *dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1) +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom ==0){value = (y1+y2)/2}else{value = num/denom}#avoid numerical errors
  return(value)
}

exp_unrevealed_1 <- function(pi,mu0,mu1,t,delta0,delta1){
  y1  = -log(t)
  y2 = -log(1-t)
  num = y1*pi*dens(t,mu1,delta1) + y2*pi*dens(1-t,mu1,delta1)
  denom =pi*dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1) +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom ==0){value = (y1+y2)/2}else{value = num/denom}#avoid numerical errors
  return(value)
}


order_prob <- function(pi,mu0,mu1,t,delta0,delta1){
  h11 = delta1*(t==1/2)+(1-delta1)*(t!=1/2)*(t^{1/mu1-1}/mu1)
  h10 = delta1*(t==1/2)+(1-delta1)*(t!=1/2)*((1-t)^{1/mu1-1}/mu1)
  h01 = delta0*(t==1/2)+(1-delta0)*(t!=1/2)*(t^{1/mu0-1}/mu0)
  h00 = delta0*(t==1/2)+(1-delta0)*(t!=1/2)*((1-t)^{1/mu0-1}/mu0)
  num = pi*h11
  denom = pi*h11+pi*h10+(1-pi)*(h00+h01)
  if(denom ==0){value =0.5}else{value = num/denom}#avoid numerical errors
  return(value)
}
