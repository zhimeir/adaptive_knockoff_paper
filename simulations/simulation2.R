#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}

####################################
## Libraries and sources
####################################
library("adaptMT")
library("splines")
library("knockoff")
library("SNPknock")
library("dplyr")
library("corpcor")
library("glmnet")
library("MASS")
library("R.matlab")
library("gam")
library("SNPknock")
library("randomForest")
library("devtools")
library("RCurl")
library("ggplot2")
library("httr")
library("mgcv")
library("flare")
library("expm")
library("rdist")
source("https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R")
source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')

source_gitfile <- function(filename){
  url = paste0("https://raw.githubusercontent.com/zhimeir/adaptive_knockoff_paper/master/",filename,".R")
  script = GET(url = url)
  script = content(script,"text")
  eval(parse(text = script),envir= .GlobalEnv)
}

file_vec <- c("utils/all_other_methods","R/adaptive_knockoff","R/filter_EM","R/filter_gam","R/filter_glm","R/filter_randomForest")
getfile <- sapply(file_vec,source_gitfile)


####################################
## Parameters
####################################
set.seed(24601)
n = 1000
p = 1600 # 30-by-30

alphalist = seq(0.3,0.01,-0.01)
Sigma = diag(rep(1,p))
amp = 25
settingName = "./results/simulation2"

####################################
## determining the signals
####################################
nonzero = c()
cor_hypothesis = expand.grid(-20:19,-20:19)
for (i in 1:p){
  xh = (cor_hypothesis[i,1])/9
  yh =(cor_hypothesis[i,2])/9
  dist = (xh^2+yh^2-2)^4-xh^3*yh^5
  dist3 = ((xh+1)^2+(yh-10/8)^2-0.015)
  dist4 = ((xh-1)^2+(yh+10/8)^2-0.015)
  
  if(dist<0 ) nonzero = c(nonzero,i)
  if(dist3<0 ) nonzero = c(nonzero,i)
  if(dist4<0 ) nonzero = c(nonzero,i)
}

# covariance matrix
Sigma = exp(-3*pdist(cor_hypothesis, metric = "euclidean", p = 2)^2)
k = length(nonzero)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
y.sample = function(X) rbinom(n,1,exp(X %*% beta0)/(1+exp(X %*% beta0)))

####################################
## Generating data
####################################
set.seed(ParamsRowIndex)
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
y = y.sample(X)

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
res = filter_glm(W,z,alpha =alphalist,offset=1)
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
res = filter_gam(W,z,alpha =alphalist,offset=1,df=5,reveal_prop = 0.2)
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
res = filter_randomForest(W,z,alpha =alphalist,offset=1)
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
res = filter_EM(W,z,alpha =alphalist,offset=1,s0=0.01,cutoff = 0)
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


