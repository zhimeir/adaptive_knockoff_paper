####################################
## Parameters
####################################
set.seed(24601)
n = 1000
p = 900 
k = 150

alphalist = seq(0.3,0.01,-0.01)
rho = 0.5

Sigma = toeplitz(rho^(0:(p-1)))
amp = 3.5
sigprob = rep(0,p)
sigprob[1:300] = 1/(1:300)^2/(sum(1/(1:300)^2))
nonzero = sample(1:p,k,prob = sigprob)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
y.sample = function(X) X%*%beta0+rnorm(n,0,1)
#settingName = paste0(".//scratch/users/zren/adaptiveKnockoff/results/",settingName)
settingName = paste0("./results/",settingName)

####################################
## Simulations
####################################
K=5;M=3;  # Number of possible states for each variable
# Marginal distribution for the first variable
pInit = rep(1/K,K)
# Create p-1 transition matrices
Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
pEmit = array(stats::runif(p*M*K),c(p,M,K))
for(j in 1:p) { pEmit[j,,] = pEmit[j,,] / rowSums(pEmit[j,,]) }



####################################
## Generating data
####################################
set.seed(ParamsRowIndex)
X = sampleHMM(pInit, Q, pEmit, n=n)
y = y.sample(X)

####################################
## Computing p values
####################################
#mdl = glm(y~X,family = binomial)
mdl = lm(y~X)
#mdl = SSLasso(X,y)
#pvals = mdl$pvals
pvals = summary(mdl)$coefficients[-1,4] 
#hist(pvals[-nonzero])
#hist(pvals)
#plot(-log(pvals))
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
z <-1:p
pi.formulas <- paste0("ns(z, df = ", 6:10, ")")
mu.formulas <- paste0("ns(z, df = ", 6:10, ")")
res <- adapt_gam(x = data.frame(z), pvals = pvals, pi_formulas = pi.formulas, mu_formulas = mu.formulas,alpha = alphalist)

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
Xk = knockoffHMM(X, pInit, Q,pEmit,seed = ParamsRowIndex+24601)
mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
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
res = adaptiveKnockoff::filter_glm(W,z,alpha =alphalist,offset=1,reveal_prop = 0.8)
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
res = adaptiveKnockoff::filter_gam(W,z,alpha =alphalist,offset=1,reveal_prop = 0.8)
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


savedir = paste(settingName,'/adakn_gam',as.character(ParamsRowIndex),'.mat',sep = "")
writeMat(savedir,fdp = fdp,power = power,res = res)



####################################
## Adaknockoff with random forest
####################################
res = adaptiveKnockoff::filter_randomForest(W,z,alpha =alphalist,offset=1,reveal_prop = 0.8)
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
#res = adaknockoff_EM(W,z,alpha =alphalist,offset=1,df = 3)
res =adaptiveKnockoff::filter_EM(W,z,alpha =alphalist,offset=1,df = 2)
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


