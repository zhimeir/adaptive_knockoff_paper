calculate.fdphat = function(W,tau,offset = 1){
  p = length(W)
  fdphat = (offset+sum(W<=-tau))/sum(W>=tau)
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
  if(denom == 0 ){value =0.5}else{value = num/denom}
  return(value)
}

prob_unrevealed <- function(pi,mu0,mu1,t,delta0,delta1){
  num = pi*dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1)
  denom = num +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if((denom) ==0){value = 0.5}else{value = num/denom}
  return(value)
}


exp_unrevealed_1 <- function(pi,mu0,mu1,t,delta0,delta1){
  y1  = -log(t)
  y2 = -log(1-t)
  num = y1 * pi *dens(t,mu1,delta1) + y2*pi*dens(1-t,mu1,delta1)
  denom =pi *dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1) +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom ==0){value = (y1+y2)/2}else{value = num/denom}
  return(value)
}

exp_unrevealed_0 <- function(pi,mu0,mu1,t,delta0,delta1){
  y1  = -log(t)
  y2 = -log(1-t)
  num = y1 *( 1-pi) *dens(t,mu0,delta0) + y2*(1-pi)*dens(1-t,mu0,delta0)
  denom =pi *dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1) +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom ==0){value = (y1+y2)/2}else{value = num/denom}
  return(value)
}


order_prob <- function(pi,mu0,mu1,t,delta0,delta1){
  h11 = delta1*(t==1/2)+(1-delta1)*(t!=1/2)*(t^{1/mu1-1}/mu1)
  h10 = delta1*(t==1/2)+(1-delta1)*(t!=1/2)*((1-t)^{1/mu1-1}/mu1)
  h01 = delta0*(t==1/2)+(1-delta0)*(t!=1/2)*(t^{1/mu0-1}/mu0)
  h00 = delta0*(t==1/2)+(1-delta0)*(t!=1/2)*((1-t)^{1/mu0-1}/mu0)
  num = pi*h11
  denom = pi*h11+pi*h10 +(1-pi)*(h00+h01)
  if(denom ==0){value =0.5}else{value = num/denom}
  return(value)
}

adaknockoff_EM <- function(W,z,df = 3,alpha = 0.1,offset = 1,s0 = 5e-3,mute = TRUE){
  p = length(W)
  z = as.matrix(z)
  all_id = 1:p
  R = 1

  # Initializing the output
  rejs = vector("list",length(alpha))
  nrejs =  rep(0,length(alpha))
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()


  # Algorithm initialization
  W_abs = abs(W)
  W_revealed = W_abs


  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)

  revealed_id = which(W_abs<=tau)


  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  W_revealed[revealed_id] = W[revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }

  pi = rep(sum(W>0)/p,p)
  delta0 = sum(W==0)/p*(1-mean(pi))
  delta1 = sum(W==0)/p*(mean(pi))
  t = logis(W_revealed)
  mu_0 =  rep(-log(logis(mean(W_abs[W<0]))),p)
  mu_1 = rep(-log(logis(mean(W_abs[W>0]))),p)
  mu_1[W==0] = log(2)
  mu_0[W==0] = log(2)  
  
  H = rep(1e-10,p)
  y0 = -log(t)
  y1 = -log(t)
  count = 0

  for (talpha in 1:length(alpha)){

    fdr = ordered_alpha[talpha]
    for (i in 1:length(unrevealed_id)){
      # EM algorithm
      for (r in 1:R){
        # E step
          # revealed part
          H[revealed_id] = sapply(1:length(revealed_id),function(j) prob_revealed(pi[revealed_id[j]],mu_0[revealed_id[j]],mu_1[revealed_id[j]],t[revealed_id[j]],delta0,delta1))
          y1[revealed_id] = -log(t[revealed_id])
          y0[revealed_id] = -log(t[revealed_id])
          
          # unrevealed part
          H[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) prob_unrevealed(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))
          H = pmax(H,1e-10)
          H = pmin(H,1-1e-10)
          y1[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) exp_unrevealed_1(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))/H[unrevealed_id]
          y0[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) exp_unrevealed_0(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))/(1-H[unrevealed_id])
          

          #if(count == 0){
        #M step
          
          if(length(unrevealed_id)<200){
            if(dim(z)[2] ==1){
              mdl = gam(log(H/(1-H))~ns(z,df))
              pi = logis(mdl$fitted.values)  
            }else{
            mdl = gam(log(H/(1-H))~s(z[,1],z[,2]))
            pi = logis(mdl$fitted.values)
            }
            # mdl =gam(H~s(z[,1],z[,2]),family = betar(link = "logit"))
            # pi =mdl$fitted.values
            if(dim(z)[2]==1){mdl =gam(y0[t!=1/2]~ns(z[t!=1/2],df),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
            }else{
              mdl =gam(y0[t!=1/2]~s(z[t!=1/2,1],z[t!=1/2,2]),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
            }
            mu_0[t!=1/2] =mdl$fitted.values
            
            if(dim(z)[2]==1){mdl =gam(y1[t!=1/2]~ns(z[t!=1/2],df),weights = (H[t!=1/2]),family = Gamma(link = "log"))
            }else{
              mdl =gam(y1[t!=1/2]~s(z[t!=1/2,1],z[t!=1/2,2]),weights = (H[t!=1/2]),family = Gamma(link = "log"))
            }
            mu_1[t!=1/2] =mdl$fitted.values  
          }else{
            mdl = randomForest(y= H,x=as.matrix(z))
            pi = mdl$predicted  
            mdl = randomForest(y = y0[t!=1/2],x= as.matrix(z[t!=1/2,]))
            mu_0[t!=1/2] =mdl$predicted
            mdl = randomForest(y = y1[t!=1/2],x= as.matrix(z[t!=1/2,]))
            mu_1[t!=1/2] =mdl$predicted
          }
            
          # if(dim(z)[2]==1){mdl =gam(H~ns(z,10),family = Gamma(link = "logit"))
          # }else{
          #   mdl =gam(H~s(z[,1],z[,2]),family = betar(link = "logit"))
          # }
          # pi =mdl$fitted.values
          
          
          
          
          
          
          #}
          delta0 = sum((1-H)*(t==1/2))/(sum((1-H)*(t==1/2))+sum((1-H)*(t!=1/2)))
          delta1 = sum((H)*(t==1/2))/(sum((H)*(t==1/2))+sum((H)*(t!=1/2)))
          count = count+1
          if(count>=20) {count  = 0}
      }


      horder = sapply(1:length(unrevealed_id),function(j) order_prob(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))
      ploth = rep(0,p)
      ploth[unrevealed_id] = sign(W[unrevealed_id])*horder
      plot(ploth)
      
      #ind.min = which(H[unrevealed_id] == min(H[unrevealed_id]))
      ind.min = which(horder == min(horder))
      if(length(ind.min)==1){
        ind.reveal = ind.min
      }else{
        ind.reveal = ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal = unrevealed_id[ind.reveal]
      revealed_id = c(revealed_id,ind.reveal)
      unrevealed_id = all_id[-revealed_id]
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      if(mute == FALSE) print(fdphat)
      if(fdphat<=fdr){break}
      if(fdphat == Inf){break}
      W_revealed[ind.reveal] = W[ind.reveal]
      t = logis(W_revealed)


    }
    if(fdphat == Inf){
      break}
    if(fdphat<=fdr){
      rej = which(W>=tau)
      rejs[[talpha]] = rej
      nrejs[talpha] = length(rej)
    }
  }

  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id)
  return(result)


}

