################################################################################
####  THIS FILE DEMONSTRATES THE TWO-STAGE ALGORITHM PROPOSED IN THE PAPER  #### 
####  "BAYESIAN VARIABLE SELECTION FOR ULTRAHIGH-DIMENSIONAL SPARSE LINEAR  ####
#### MODEL" WITH AN APPLICATION.                                            ####
################################################################################

rm(list = ls())
library(Matrix)
library(mvtnorm)
library(elasticnet)
library(glmnet)
library(SIS)
library(ncvreg)
library(BLR)

source("functions.R")
################################ Fixed parameters ##############################
n = 100 
p = 1000 
DF = 0 
rept = 1
tm = ifelse(p <= 500, 5, 10)
u = rep(2, tm)

beta = rep(0, p)
beta[1] = 1
beta[(2:(tm+1))] = u
permute = sample(c(2:p),(p-1),replace=F)
beta[permute] = beta[-1]
beta_star = which(beta != 0)
indx_comp = which(beta == 0)
beta_val = beta[beta_star]
count_stp2 = count_stp1 = rep(0,rept)

################## construction of sigma matrix ################################
########### AR structure
r = 0.5
sig = matrix(0,p,p)
for (i in 1:p){
  for (j in i:p){
    if (j == i) sig[i,j] = 1
    if (j > i) sig[i,j] = r^(abs(j-i))
    sig[j,i] = sig[i,j]
  }
}
Sigma = sig

################################################################################
############################# Constructing the X matrix ########################
################################################################################
x = matrix(data = 1, nrow = n, ncol = p)
x[(1:n),] = rmvnorm(n, mean = array(data = 0, dim = p), sigma = Sigma)
  #x[(1:n),]=rmvt(n,delta=rep(0,p),sigma=Sigma,df=3)
x = center(x)
x[,1] = 1
MU_0 = x %*% beta;
  #DF=2  #0,2
if(DF == 0) {en = rmvnorm(1, mean = array(data = 0,dim = n), sigma = diag(n))}
if(DF == 2) {en = rmvt(1, delta = rep(0, n), sigma = diag(n), df = DF)}
  
############################### Construction of y ##############################
Y = y = MU_0 + t(en)
X = x
  #  if(r == 1) print("data done")
  
############################################################################
###############################    METHOD   ################################
############################################################################
g = n * p
N = 10
D = floor(n/4) #25, floor(n/2)  #subjective
start = Sys.time()
gmm_super = step1mcmc(X, y, D, g, N)
end = Sys.time()
runtime.step1 = difftime(end, start, unit = "secs")
S = which(gmm_super == 1)
beta_n = beta[S]
x_n = X[,S]
df_stp1 = length(setdiff(beta_star,S))
if(df_stp1 == 0)
{
  count_stp1[r] = 1
}
gmm_star = rep(1, length(beta_n))
gmm_star[which(beta_n == 0)] = 0
gmm_star = as.vector(gmm_star)

#source("functions4step2.r")
g = D^2 
Repl = 6000
w = 2
start = Sys.time()
result = MCMC(S, x_n, gmm_star, y, df_stp1, g, w, Repl)  
end = Sys.time()
runtime.step2 = difftime(end, start, unit = "secs")
#posterior[r]=result[1]
Probability = result[[1]]
mode.model = result[[2]]
DIST = (length(setdiff(S[mode.model == 1], beta_star)) + 
          length(setdiff(beta_star, S[mode.model == 1])))
if(DIST == 0)
{
  count_stp2[r] = 1
}


print(date())
cat("---------------------------------------------------","\n")
cat("(n, p) =", c(n, p),"\n")
cat("The indices of the TRUE covariate is =", beta_star,"\n")
cat("The indices of the first-stage-selected model is =", S,"\n")
cat("The number of true covariates excluded after screening =", df_stp1,"\n")
cat("The time (in seconds) taken to calculate the first-stage is =", runtime.step1,"\n")

cat("The indices of the MODE model is =", S[mode.model == 1],"\n")
cat("The empirical probability of the MODE model is =", Probability,"\n")
cat("The time (in seconds) taken to calculate the second-stage is =", runtime.step2,"\n")
cat("The Hamming distance of TRUE model and MODE model is =", DIST,"\n")
cat("---------------------------------------------------","\n")
