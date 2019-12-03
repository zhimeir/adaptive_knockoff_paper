####################################
## Parameters
####################################
set.seed(24601)
n = 1000
p = 1600 # 30-by-30
#k=100 # 6-by-6

alphalist = seq(0.3,0.01,-0.01)

#rho = 0.1

#Sigma = toeplitz(rho^(0:(p-1)))
Sigma = diag(rep(1,p))
amp = 25
#settingName = paste0("/scratch/users/zren/adaptiveKnockoff/results/",settingName)
settingName = paste0("./results/",settingName)
####################################
## Ploting the signals
####################################
nonzero = c()
cor_hypothesis = expand.grid(-20:19,-20:19)
for (i in 1:p){
  xh = (cor_hypothesis[i,1])/9
  yh =(cor_hypothesis[i,2])/9
  dist = (xh^2+yh^2-2)^4-xh^3*yh^5
  #dist2 = (xh^2+yh^2-0.05)
  dist3 = ((xh+1)^2+(yh-10/8)^2-0.015)
  dist4 = ((xh-1)^2+(yh+10/8)^2-0.015)
  
  if(dist<0 ) nonzero = c(nonzero,i)
  #if(dist2<0 ) nonzero = c(nonzero,i)
  if(dist3<0 ) nonzero = c(nonzero,i)
  if(dist4<0 ) nonzero = c(nonzero,i)
}

# and the covaiance matrix should have corresponding structure
# for (i in 2:p){ 
#   for(j in 1:(i-1)){
#     #Sigma[i,j] = rho*(abs(cor_hypothesis[i,1]-cor_hypothesis[j,1])+abs(cor_hypothesis[i,2]-cor_hypothesis[j,2])==1)
#     Sigma[i,j] = exp(-3*sum((cor_hypothesis[i,] - cor_hypothesis[j,])^2))
#     Sigma[j,i] = Sigma[i,j]
#   }
# }
Sigma = exp(-3*pdist(cor_hypothesis, metric = "euclidean", p = 2)^2)

#group = rep(1,p)
#group[nonzero]= 19
#color = rep("black",p)
#color[nonzero] = "red"

#plot(x = cor_hypothesis[,1],y= cor_hypothesis[,2],pch = group,col = color)
k = length(nonzero)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
y.sample = function(X) rbinom(n,1,exp(X %*% beta0)/(1+exp(X %*% beta0)))
#y.sample = function(X) X%*%beta0+rnorm(n,0,1)


####################################
## Simulations
####################################

####################################
## Generating data
####################################
set.seed(ParamsRowIndex)
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
y = y.sample(X)

####################################
## Computing p values
####################################
#mdl = glm(y~X[,1],family = binomial)
 #mdl = lm(y~X)
#pvals = summary(mdl)$coefficients[-1,4]
#hist(pvals[-nonzero])
get_pval = function(x){
  mdl = glm(y~x,family = binomial)
  pval = summary(mdl)$coefficients[-1,4]
  return(pval)
}
get_pval(X[,1])
pvals <- apply(X,MARGIN = 2,FUN = get_pval)

####################################
## BHq
####################################
fdp = c()
power = c()
for(i in 1:length(alphalist)){
  rej = BH_method(pvals, alphalist[i])
  fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
  power = c(power,sum(beta0[rej]!=0)/max(k,1))
}
savedir = paste(settingName,'/BH',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir, power = power, fdp = fdp,nonzero = nonzero)

####################################
## Storey method
####################################
thr = 0.5 #(as in Rina and Li)
fdp = c()
power = c()
for(i in 1:length(alphalist)){
  rej = Storey_method(pvals, alphalist[i], thr)
  fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
  power = c(power,sum(beta0[rej]!=0)/max(k,1))
}
savedir = paste(settingName,'/storeyBH',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir, power = power, fdp = fdp,nonzero = nonzero)

####################################
## SABHA method
####################################
tau = 0.5;eps = 0.1;#as in Rina and Li
qhat = Solve_q_step(pvals,tau,eps)
fdp = c()
power = c()

for(i in 1:length(alphalist)){
  rej =SABHA_method(pvals, qhat, alphalist[i], tau)
  fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
  power = c(power,sum(beta0[rej]!=0)/max(k,1))
}
savedir = paste(settingName,'/SABHA',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir, power = power, fdp = fdp,nonzero = nonzero)


####################################
## Adaptive Seqstep
####################################
thr2= 0.5
fdp = c()
power = c()

for(i in 1:length(alphalist)){
  thr1 = alphalist[i]
  rej =Adaptive_SeqStep_method(pvals, alphalist[i], thr1, thr2)
  fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
  power = c(power,sum(beta0[rej]!=0)/max(k,1))
}
savedir = paste(settingName,'/seqstep',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir, power = power, fdp = fdp,nonzero = nonzero)

####################################
## AdaPT
####################################
z <- data.frame(xaxis = cor_hypothesis[,1],yaxis = cor_hypothesis[,2])
formulas <- c(paste0("ns(xaxis, df = ", 6:10, ")"),paste0("ns(yaxis, df = ", 6:10, ")"))
pi_formula <- mu_formula <- "s(xaxis, yaxis)"

res <- adapt_gam(x = z, pvals = pvals, pi_formulas = pi_formula, mu_formulas = mu_formula,alpha = alphalist)


fdp = c()
power = c()
for (i in 1:length(alphalist)){
  if(res$nrejs[i]>0){
    rej =res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))
  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}
savedir = paste(settingName,'/adapt',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir, power = power, fdp = fdp,nonzero = nonzero)

####################################
## Vanilla  knockoff
####################################
Xk = create.gaussian(X,rep(0,p),Sigma)
mdl = cv.glmnet(cbind(X,Xk),y,alpha=1,family = "binomial")
cvlambda = mdl$lambda.min
beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
T = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T,T_tilde)
W = T-T_tilde

fdp = c()
power = c()
for (i in alphalist ){
  tau = knockoff.threshold(W,fdr = i,offset = 1)
  rej = which(W>=tau)
  fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
  power = c(power,sum(beta0[rej]!=0)/max(k,1))
}
savedir = paste(settingName,'/vanilla',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir, power = power, fdp = fdp,W = W,nonzero = nonzero)

####################################
## Adaknockoff with glm
####################################
z = as.matrix(cor_hypothesis)
res = adaptiveKnockoff::filter_glm(W,z,alpha =alphalist,offset=1)
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))

  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}


savedir = paste(settingName,'/adakn_glm',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir,fdp = fdp,power = power,res = res)


####################################
## Adaknockoff with gam
####################################
# res = adaptiveKnockoff::filter_gam(W,z,alpha =alphalist,offset=1, df = 2)
# fdp = c()
# power = c()
# for (i in 1:length(alphalist)){
#   if(res$nrejs[[i]]>0){
#     rej  = res$rejs[[i]]
#     fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
#     power = c(power,sum(beta0[rej]!=0)/max(k,1))
# 
#   }else{
#     fdp = c(fdp,0)
#     power = c(power,0)
#   }
# }
# 

savedir = paste(settingName,'/adakn_gam',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir,fdp = fdp,power = power,res = res)



####################################
## Adaknockoff with random forest
####################################
res = adaptiveKnockoff::filter_randomForest(W,z,alpha =alphalist,offset=1)
fdp = c()
power = c()
for (i in 1:length(alphalist)){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))

  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}


savedir = paste(settingName,'/adakn_rf',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir,fdp = fdp,power = power,res = res)

####################################
## Adaknockoff with EM algorithm
####################################
#res = adaknockoff_EM(W,z,alpha =alphalist,offset=1,s0=0.01,mute = FALSE)
res = adaptiveKnockoff::filter_EM(W,z,alpha =alphalist,offset=1,s0=0.01,df=2,cutoff = 100)
fdp = c()
power = c()
for (i in 1:length(alphalist )){
  if(res$nrejs[[i]]>0){
    rej  = res$rejs[[i]]
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))

  }else{
    fdp = c(fdp,0)
    power = c(power,0)
  }
}


savedir = paste(settingName,'/adakn_em',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir,fdp = fdp,power = power)


