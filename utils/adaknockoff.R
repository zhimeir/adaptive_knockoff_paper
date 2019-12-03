calculate.fdphat = function(W,tau,offset = 1){
  p = length(W)
  fdphat = (offset+sum(W<=-tau))/sum(W>=tau)
  return(fdphat)
}



##############################
###    glmnet version      ###
##############################

adaknockoff_glmnet <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  #nnz = which(W!=0)
  #W = W[nnz]
  #z= z[nnz]
  p = length(W)
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  W_abs = abs(W)
  W_sign = as.numeric(W>0)
  revealed_sign = rep(1,p)
  tau = rep(s0,p)
  all_id = 1:p
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]}else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl = cv.glmnet(cbind(z,abs(W)),revealed_sign,family = "binomial")
      fitted.pval = predict(mdl,cbind(z,abs(W)),s= "lambda.min",type = "response")
      fitted.pval = fitted.pval[unrevealed_id]
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
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = rej.path,unrevealed_id = unrevealed_id,selected.pval)
  return(result) 
}

##############################
### random forest version  ###
##############################


adaknockoff_rf <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  p = length(W)
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  tau.sel = c()
  acc = c()
  
  W_abs = abs(W)
  W_sign = sign(W)/2+1/2
  revealed_sign = rep(1,p)
  tau = rep(s0,p)
  all_id = 1:p
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id

  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
  
    for (i in 1:length(unrevealed_id)){
      mdl = randomForest(y = as.factor(revealed_sign),x = cbind(W_abs,z),norm.votes = TRUE)
      fitted.pval = mdl$votes[,ncol(mdl$votes)]
      fitted.pval = fitted.pval[unrevealed_id]
      predicted.sign = fitted.pval>0.5
      acc = c(acc, sum(predicted.sign == W_sign[unrevealed_id])/length(unrevealed_id))
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
    tau.sel = c(tau.sel, length(revealed_id))
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,sel.prob = selected.pval,tau = tau.sel,acc = acc)
  return(result) 
}


##############################
# linear regression version  #
##############################

adaknockoff_lm <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  
  p = length(W)
  all_id = 1:p
  
  # Initializing the output
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  tau.sel = c()
  acc = c()
  
  
  # Algorithm initialization
  W_abs = abs(W)
  W_sign = sign(W)/2+1/2
  revealed_sign = rep(1,p)
  
  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)
  
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl = glm(revealed_sign~cbind(z,abs(W)),family = "binomial")
      fitted.pval = predict(mdl,type = "response")
      fitted.pval = fitted.pval[unrevealed_id]
      predicted.sign = fitted.pval>0.5
      acc = c(acc, sum(predicted.sign == W_sign[unrevealed_id])/length(unrevealed_id))
      
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
    tau.sel = c(tau.sel, length(revealed_id))
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,sel.prob = selected.pval,tau = tau.sel,acc = acc)
  return(result) 
}


##############################
# gam regression version  ####
##############################

adaknockoff_gam <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
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
  tau = rep(s0,p)
  all_id = 1:p
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl  =gam(revealed_sign~cbind(z,abs(W)),family = "binomial")
      fitted.pval = mdl$fitted.values
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


##############################
# linear regression version  #
##############################

adaknockoff_lm_beta <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  
  p = length(W)
  all_id = 1:p
  
  # Initializing the output
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  pos_prob = c()
  
  # Algorithm initialization
  W_abs = abs(W)
  W_sign = sign(W)/2+1/2
  revealed_sign = rep(1,p)
  
  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)
  
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  neg_num = sum(W<=-tau)
  pos_num = sum(W>=tau)
  revealed_sign[-revealed_id] = pos_num/(pos_num+neg_num)
  pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl = glm(revealed_sign~cbind(z,abs(W)),family = "binomial")
      fitted.pval = predict(mdl,type = "response")
      fitted.pval = fitted.pval[unrevealed_id]
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
      pos_num = sum(W[unrevealed_id]>0)
      neg_num = sum(W[unrevealed_id]<0)
      pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
      revealed_sign[unrevealed_id] = pos_num/(pos_num+neg_num)
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      print(fdphat)
      if(fdphat<=fdr){break}
    }
    
    rej = which(W>=tau)  
    rejs[[talpha]] = rej
    nrejs[[talpha]] = length(rej)
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,selected.pval,numremain = length(unrevealed_id),pos_prob = pos_prob)
  return(result) 
}

##############################
# gam version  ###############
##############################

adaknockoff_gam_beta <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  
  p = length(W)
  all_id = 1:p
  
  # Initializing the output
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  pos_prob = c()
  
  # Algorithm initialization
  W_abs = abs(W)
  W_sign = sign(W)/2+1/2
  revealed_sign = rep(1,p)
  
  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)
  
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  neg_num = sum(W<=-tau)
  pos_num = sum(W>=tau)
  revealed_sign[-revealed_id] = pos_num/(pos_num+neg_num)
  pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl = gam(revealed_sign~cbind(z,abs(W)),family = "binomial")
      fitted.pval = mdl$fitted.values
      fitted.pval = fitted.pval[unrevealed_id]
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
      pos_num = sum(W[unrevealed_id]>0)
      neg_num = sum(W[unrevealed_id]<0)
      pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
      revealed_sign[unrevealed_id] = pos_num/(pos_num+neg_num)
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      print(fdphat)
      if(fdphat<=fdr){break}
    }
    
    rej = which(W>=tau)  
    rejs[[talpha]] = rej
    nrejs[[talpha]] = length(rej)
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,selected.pval,numremain = length(unrevealed_id),pos_prob = pos_prob)
  return(result) 
}



##############################
# random forest version  #####
##############################

adaknockoff_rf_beta <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  
  p = length(W)
  all_id = 1:p
  
  # Initializing the output
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  pos_prob = c()
  tau.sel = c()
  acc = c()
  
  # Algorithm initialization
  W_abs = abs(W)
  W_sign = sign(W)/2+1/2
  revealed_sign = rep(1,p)
  
  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)
  
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  neg_num = sum(W<=-tau)
  pos_num = sum(W>=tau)
  revealed_sign[-revealed_id] = pos_num/(pos_num+neg_num)
  pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      mdl = randomForest(y = as.factor(revealed_sign),x = cbind(W_abs,z),norm.votes = TRUE)
      dm = ncol(mdl$votes)
      fitted.pval = mdl$votes[,dm]+ mdl$votes[,dm-1]
      fitted.pval = fitted.pval[unrevealed_id]
      predicted.sign = fitted.pval>0.5
      acc = c(acc,sum(predicted.sign == W_sign[unrevealed_id])/length(unrevealed_id))
      
      
      
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
      pos_num = sum(W[unrevealed_id]>0)
      neg_num = sum(W[unrevealed_id]<0)
      pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
      revealed_sign[unrevealed_id] = pos_num/(pos_num+neg_num)
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
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,sel.prob = selected.pval,tau = tau.sel,pos_prob = pos_prob,
            acc = acc)
  return(result) 
}

##############################
########## CNN version  #####
##############################

adaknockoff_cnn <- function(W,z,alpha =0.1,offset=1,s0 = 1e-3){
  
  p = length(W)
  all_id = 1:p
  
  # Initializing the output
  rejs = list()
  nrejs = list()
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  selected.pval = c()
  pos_prob = c()
  
  # Algorithm initialization
  W_abs = abs(W)
  W_sign = as.numeric(W>0)
  revealed_sign = rep(1,p)
  
  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)
  
  revealed_id = which(W_abs<=tau)
  revealed_sign[revealed_id] = W_sign[revealed_id]
  neg_num = sum(W<=-tau)
  pos_num = sum(W>=tau)
  #revealed_sign[-revealed_id] = pos_num/(pos_num+neg_num)
  pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
  
  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  rej.path = c(rej.path,revealed_id)
  }else{
    unrevealed_id = all_id
  }
  
  
  for (talpha in 1:length(alpha)){
    
    fdr = ordered_alpha[talpha]
    
    for (i in 1:length(unrevealed_id)){
      x_train <- cbind(W_abs[revealed_id],z[revealed_id])
      y_train <- revealed_sign[revealed_id]
      y_train <- to_categorical(y_train,2)
      
      x_test <-  cbind(W_abs[unrevealed_id],z[unrevealed_id])
      y_test <- W_sign[unrevealed_id]
      y_test <- to_categorical(y_test,2)
      
      model <- keras_model_sequential() 
      model %>% 
        layer_dense(units = 10, activation = 'relu', input_shape = c(2)) %>% 
        layer_dense(units = 10, activation = 'relu') %>% 
        layer_dense(units = 2, activation = 'softmax')
      
      
      model %>% compile(
        loss = 'categorical_crossentropy',
        optimizer = optimizer_adam(),
        metrics = c('accuracy')
      )
      
      history <- model %>% fit(
        x_train, y_train, 
        epochs = 30, batch_size = 50, 
        validation_split = 0,
        verbose = 0
      )
      
      output = model %>% predict(x_test)
      #print(model%>%evaluate(x_test,y_test))
      fitted.pval = output[,2]
      #fitted.pval = fitted.pval[unrevealed_id]
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
      pos_num = sum(W[unrevealed_id]>0)
      neg_num = sum(W[unrevealed_id]<0)
      pos_prob = c(pos_prob,pos_num/(pos_num+neg_num))
      #revealed_sign[unrevealed_id] = pos_num/(pos_num+neg_num)
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      print(fdphat)
      if(fdphat<=fdr){break}
    }
    
    rej = which(W>=tau)  
    rejs[[talpha]] = rej
    nrejs[[talpha]] = length(rej)
  }
  
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,selected.pval,numremain = length(unrevealed_id),pos_prob = pos_prob)
  return(result) 
}



