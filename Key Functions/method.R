library(HIMA)
library(MultiMed)
library(MediationFDR)
source("mymethod.R")

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
