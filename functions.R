################################################################################
############################### FUNCTIONS ######################################
################################################################################
#################
### log gamma ###
#################
#' @description calculates the logarithm of gamma function of a postive number s
#' @param s
#' @returns log(Gamma(s)) 
log_gamma=function(s)
{
  R = (log(2*pi) - log(s))/2 + s * (log(s + 1/(12*s - 1/(10 * s))) - 1)
  return(R)
}
###################
### standardize ###
###################
#' @description standardizes each column of a matrix 
#' @param X a nXp matrix whose columns are standardized
#' @returns The nXp matrix with standardized entries
center = function(X)
{
  X_star = apply(X, 2, function(u) {return((u - mean(u))/sd(u))})
  return(X_star)
}
######################################
### log marginal: sigma known case ###
######################################
#' @description calculates the marginal likelihood for the model gamma in sigma 
#'              known case 
#' @param gmm index set corresponding to the model gamma
#' @param x_star the standardized full design matrix
#' @param y response vector
#' @param g choice of g
#' @returns the logarithm of the marginal likelihhod of the gamma-th model
lg_marg = function(gmm, x_star, y, g)
{
  pp = length(which(gmm == 1))
  n = length(y)
  p = length(gmm)
  sigma2 = sum((y - mean(y))^2)/n
  x_find = x_star[ , which(gmm == 1)]
  x_p = as.matrix(x_find[,-1]);
  if(pp == 1)
  {
    lmarg = 0
  }
  if(pp>1)
  {
    Pk = (t(y)) %*% (x_p %*% (solve((diag(pp - 1)/g) + 
                                      t(x_p) %*% x_p)) %*% (t(x_p))) %*% y
    lmarg = - (pp-1) * (log(p - 1) + log(g)/2) - log(det(diag(pp - 1)/g + 
                                          (t(x_p) %*% x_p)))/2 + Pk/(2 * sigma2)
  }
  return(lmarg)
}
########################################
### log marginal: sigma unknown case ###
########################################
#' @description calculates the marginal likelihood for the model gamma in sigma 
#'              known case 
#' @param gmm index set corresponding to the model gamma
#' @param x_star the standardized full design matrix
#' @param y response vector
#' @param g choice of g
#' @returns the logarithm of the marginal likelihhod of the gamma-th model
lg_su_marg=function(gmm, x_star, y, g)
{
  pp = length(which(gmm == 1))
  n = length(y)
  p = length(gmm)
  # g = max((p * p),(n * n))
  sigma = 1
  x_find = x_star[ , which(gmm == 1)]
  nsy2 = t(y) %*% y - n * mean(y) * mean(y)
  if(pp == 1)
  {
    lmarg = -(n-1) * log(nsy2)/2
  }
  if(pp>1)
  {
    x_p = as.matrix(x_find[ , -1])
    Pk = (t(y)) %*% (x_p %*% (solve((diag(pp - 1)/g) + 
                                      t(x_p) %*% x_p)) %*% (t(x_p))) %*% y
    lmarg = - (pp - 1) * (log(p - 1) + log(g)/2) - log(det(diag(pp - 1)/g + 
                                (t(x_p) %*% x_p)))/2 - (n-1) * log(nsy2 - Pk)/2
  }
  return(lmarg)
}


#######################################
### FIRST STAGE SCREENING ALGORITHM ###
#######################################
#' @description This function employs the screening algorithm as provided in
#'		the first step of the proposed method
#' @param x_n the standardized design matrix
#' @param y response vector
#' @param D the choice of d_n
#' @param Repl The number of replications 
#' @param g choice of g
#' @returns a binary sequence of length p indicating the index set of 
#'	        the first stage selected model

step1mcmc=function(x_n, y, D, g, Repl)
{
  p = ncol(x_n)
  gmm_null = array(data = 0, dim = p)
  gmm_null[1] = 1
  gmm_p = gmm_null
  gmm_p[c(2:(D+1))] = 1
  # Updating each variable gamma_i
  for(i in 1:Repl)
  {
    set1 = which(gmm_p == 0)
    set2 = which(gmm_p[2:p] == 1)
    for(j in 2:(D+1))
    {
      marg_E = NULL
      marg_E[1] = lg_su_marg(gmm_p, x_n, y, g)
      gmm_prv = gmm_p
      l1 = length(set1)
      gmm_p[set2[j]] = 0
      for(k in 1:l1)
      {
        gmm_E = gmm_p
        gmm_E[set1[k]] = 1
        marg_E[k+1] = lg_su_marg(gmm_E, x_n, y, g)
      }
      slctd = which.max(marg_E)
      if(slctd == 1)
      {
        gmm_p = gmm_prv
      }
      else
      {
        gmm_p[set1[slctd-1]] = 1
      }
      if(slctd != 1)
      { set1 = set1[-(slctd-1)] }
    }
  }
  return(gmm_p)
}


#######################################
### SECOND STAGE VARIABLE SELECTION ###
#######################################
#' @description This function employs the RJMCMC algorithm as provided in
#'		the second step of the proposed method
#' @param x_n the standardized design matrix
#' @param gmm_star the index set of the true covariates within the first stage
#'           		   selected model
#' @param y response vector
#' @param df_stp1 the number of true covariates not included in the first stage
#           		  selected model
#' @param Repl The number of iterations 
#' @param g choice of g
#' @param w choice of w in the model prior
#' @returns a list of two components: (1) the post-burnin proportion of times of 
#'    	    visit to the MODE model, i.e, the most visited model; (2) the index
#'    	    set of the MODE model within the screened covariates 

MCMC = function(S, x_n, gmm_star, y, df_stp1, g, w, Repl)
{
  new_p = length(S)
  gmm = c(1, rep(0, (new_p - 1))) # The Null Model
  ################## Updating ######################  
  count = 0
  index = NULL
  for(i in 1:Repl)	# MCMC replicate starts
  {
    # if(i %% 100 == 0)
    # {
    #   print(i)
    #   print(date())
    # }
    # # Finding the highest visited model
    if(i == floor(Repl/2))
    { 
      Locate = gmm 
      index[1] = 1	
    }
    if(i > floor(Repl/2))
    {	
      if(is.vector(Locate) == TRUE)
      {	
        if(sum(abs(gmm - Locate)) != 0)	
        {	
          Locate = rbind(Locate, gmm) 
          index[2] = 1	
        }
        else
        {	
          index[1] = index[1] + 1
        }
      }
      if(is.matrix(Locate) == TRUE)
      {
        ln = nrow(Locate)
        signal = 0
        for(lm in 1:ln)
        {	
          if(sum(abs(gmm - Locate[lm,])) == 0) 
          { 
            signal = 1 
            index[lm] = index[lm] + 1;
          }
        }
        if(signal == 0)
        {	
          Locate = rbind(Locate, gmm) 
          index[ln + 1] = 1;	
        }		
      }
    }
    ## Updating each variable gamma_i
    if(((sum(abs(gmm - gmm_star)) + df_stp1) == 0) && (i > floor(Repl/2)))
    {  
      count = count + 1		
    }
    for(k in 2:new_p)
    {
      gmm_n = gmm
      gmm_n[k] = (1-gmm[k])
      binom_prob = lg_su_marg(gmm, x_n, y, g) - lg_su_marg(gmm_n, x_n, y, g) + 
                                                (gmm_n[k] - gmm[k]) * log(w - 1)  
      # exp(log(marginal_new) - log(marginal_new + marginal) + (|gmm_new| - |gmm|) log (w-1)) 
      # binom_prob = (log(marg(gmm, x_n, y)) - log(marg(gmm_n, x_n, y)))
      samp = runif(1, min = 0, max = 1)
      #   samp=(samp)^(1+T)
      crit = (log((1/samp) - 1))
      if(crit > binom_prob)
      {    
        gmm = gmm_n	
      }
    }
  }
  if(is.vector(Locate) == TRUE)
  {
	best_model = Locate
  }
  if(is.matrix(Locate) == TRUE)
  {
  	hgt_idx = which.max(index)
  	best_model = Locate[hgt_idx,]
  }
  posterior = count/(Repl - floor(Repl/2))
  rtrn = list(posterior, best_model)
  return(rtrn)
}

