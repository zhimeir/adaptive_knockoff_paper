#!/usr/bin/env Rscript

## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
paramFile <- args[1]
ParamsRowIndex <- as.integer(args[2])
params <- read.csv(paramFile, header = TRUE, stringsAsFactors = FALSE)


####################################
## Libraries
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
#library("mht")
#library("reshape")
library("devtools")
library("RCurl")
library("keras")
library("ggplot2")
library("httr")
library("flare")
library("expm")
#install_github("zhimeir/akn",subdir = "adaptiveKnockoff",auth_token = "1105343890ba2690dd0ee046fe5f5be604e34ffd",force = TRUE)
library("adaptiveKnockoff")
library("mgcv")
source("https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R")
source('https://www.stat.uchicago.edu/~rina/accumulationtests/accumulation_test_functions.R')
url = paste0("https://raw.githubusercontent.com/zhimeir/akn/master/all_other_methods.R")
script = GET(url = url,authenticate("1105343890ba2690dd0ee046fe5f5be604e34ffd",""))
script = content(script,"text")
eval(parse(text = script))

url = paste0("https://raw.githubusercontent.com/zhimeir/akn/master/adaknockoff_em.R")
script = GET(url = url,authenticate("1105343890ba2690dd0ee046fe5f5be604e34ffd",""))
script = content(script,"text")
eval(parse(text = script))

url = paste0("https://raw.githubusercontent.com/zhimeir/akn/master/lasso_inference.r")
script = GET(url = url,authenticate("1105343890ba2690dd0ee046fe5f5be604e34ffd",""))
script = content(script,"text")
eval(parse(text = script))

#settingName = "2d_linear_gaussian_twoclusters"
settingName = "2d_linear_gaussian_twoclusters"
url = paste0("https://raw.githubusercontent.com/zhimeir/akn/master/simulations/",settingName,".R")
script = GET(url = url,authenticate("1105343890ba2690dd0ee046fe5f5be604e34ffd",""))
script = content(script,"text")
eval(parse(text = script))

