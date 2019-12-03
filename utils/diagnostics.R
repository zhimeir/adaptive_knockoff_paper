plot_rejction_path <- function(mdl,W,nonzero = NULL){
  rej_path = mdl$rej.path
  p = length(W)
  
  wpath = W[rej_path]
  nnz = which(wpath!=0)
  nz = sum(W == 0)
  if(is.null(nonzero) == 0){
    group = rep(1,p)
    group[nonzero] = 2
    
    grouppath = group[rej_path]
    grouppath = grouppath[nnz]
    
    df = data.frame(x = 1:length(nnz), W=wpath[nnz])
    ggplot(df,aes(x,W))+
      geom_bar(stat = "identity", aes(fill =grouppath),alpha = c(rep(1, mdl$tau[length(mdl$tau)]-nz),rep(0.5,p-mdl$tau[length(mdl$tau)])),show.legend = FALSE)+
      geom_vline(xintercept = mdl$tau-nz,color = "red", linetype = "dashed")
  }else{
    df = data.frame(x = 1:length(nnz), W=wpath[nnz])
    ggplot(df,aes(x,W))+
      geom_bar(stat = "identity",show.legend = FALSE)+
      geom_vline(xintercept = mdl$tau-nz,color = "red")
  }
}

plot_threshold <- function(mdl,W,nonzero){
  p = length(W)
  
  k = length(mdl$nrejs)
  
  for (talpha in 1:k){
    rej = res$rejs[[talpha]]
    adathres = rep(0,p)
    upp = min(W[rej])
    low = max(abs(W[-rej]))
    for (i in 1:p){
      if(i%in%rej==1 ){
        adathres[i] = upp
      }else{
        adathres[i] = low
      }
    }
    
    pp =  ggplot()+
      geom_point(aes(x =1:p,y = W),color = group) +
      geom_smooth(aes(x = 1:p, y = adathres),span = 0.3,method = "loess",se = FALSE)+
      geom_line(aes(x = 1:p, y =adathres),color = alpha("black",0.5),linetype = "dashed")
    plot(pp)
  }
}

plot_accuracy <- function(mdl,W){
  ggplot()+
    geom_point(aes(x = 1:length(mdl$acc),y = mdl$acc))
}

