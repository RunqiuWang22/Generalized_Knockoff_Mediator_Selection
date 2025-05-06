library(HIMA)
library(MultiMed)
library(MediationFDR)
source("mymethod.R")
### data without confounding
datagen <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  beta = rep(0,p); beta[pathb] = b; beta[pathi] = b*ss;
  gamma=1
  X=rnorm(n,0,1)
  M = X %o% alpha + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  Y = X*gamma + M %*% beta + rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M)
}

### data with confounding
datagen.conf <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  XV = mvtnorm::rmvnorm(n, mean = rep(0,2), sigma = matrix(c(1,0.5,0.5,1),2,2)) # sigma is covariance matrix
  X = as.matrix(XV[,1])
  V = as.matrix(XV[,2])
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  eta1 = rep(0,p); eta1[p] = 1 ### need to consider
  M = X %*% alpha + V %*% eta1 + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  
  beta = rep(0,p); beta[pathb] = b; beta[pathi] = b*ss;
  gamma = 1
  eta2 = 1
  Y = X*gamma + M %*% beta + eta2*V+rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M, V = V)
}

datagen.nonlinear <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  beta = rep(0,p); beta[pathb] = b*ss; beta[pathi] = b;
  gamma=0
  X=rnorm(n,0,1)
  M = X %o% alpha + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  M_cos <- apply(M, 2, cos)
  Y = X*gamma + M_cos %*% beta + rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M)
}


### data real with confounding
datagen.real <- function(n, p, a, b,ss,rho,std,patha,pathb,meanM,covmat){
 # meanM: mean of M
 # covmat:covariance of Y
  pathi = intersect(patha, pathb)
  XV = mvtnorm::rmvnorm(n, mean = rep(0,2), sigma = matrix(c(1,0.5,0.5,1),2,2)) # sigma is covariance matrix
  X = as.matrix(XV[,1])
  V = as.matrix(XV[,2])
  X = ifelse(X>=median(X),1,0) ### change X into binary
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  eta1 = rep(0,p); eta1[p] = 1 ### need to consider
  M = X %*% alpha + V %*% eta1 + mvtnorm::rmvnorm(n, mean = meanM, sigma = covmat) ## sigma: covariance matrix
  
  beta = rep(0,p); beta[pathb] = b; beta[pathi] = b*ss;
  gamma = 1
  eta2 = 1
  Y = X*gamma + M %*% beta + eta2*V+rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M, V = V)
}

myest = function(X,Y,M,V=NULL,method,q){
  p=dim(M)[2]
  n=dim(M)[1]
 
  
  if (method=="MCP"){
    XV=cbind(X,V)
    ###################  MCP #########################
    Mnames <- colnames(M)
    ##pEM: get exposure/mediator relationships
    ##pMY: get mediator/outcome relationship (conditional on exposure)
    
    pEM = pMY= rep(0,p)
    for (j in 1:p){
      pEM[j] = coef(summary(lm(M[,j] ~ XV)))[2,4]
      pMY[j] = coef(summary(lm(Y~M[,j]+XV)))[2,4]
    } 
    
    medTest.FDR <- medTest.SBMH(pEM, pMY, MCP.type="FDR", t1 = 0.1, t2 = 0.1)
    out = Mnames[which(medTest.FDR < q)]
    select = as.numeric(gsub("[^0-9]", "", out))
  }
  
  if (method=="HIMA"){
    ##################### HIMA ###########################
    if (!is.null(V)) {V=as.matrix(V)}
    result = classicHIMA1(X=X,Y=Y,M=M,COV.XM=V,COV.MY=V,Bonfcut=q)
    select = as.numeric(gsub("[^0-9]", "", result$Index))
  }
  
  if (method=="HIMA2"){
    #################### HIMA2 ###########################
    result = dblassoHIMA(X=X,Y=Y,M=M,COV=V,FDRcut=q)
    select = as.numeric(gsub("[^0-9]", "", result$Index))
  }
  
  if (method=="eHIMA"){
    #################### eHIMA ###########################
    result = eHIMA1(X=X,Y=Y,M=M,COV=V,FDRcut=q)
    select = as.numeric(gsub("[^0-9]", "", result$Index))
  }
  
  if (method=="MediationFDR"){
    ################# MediationFDR ###########################
    result = MediationFDR(X,Y,M,V1=V,V2=V,q1=0.1,q2=0.1)
    select = result$med_select
  }
  
  if (method=="Mymethod"){
    ############### proposed method #########################
    select = mymethod(X=X,Y=Y,M=M,V=V,q=q,type="non-linear")
  }
  
  return(select)
}
