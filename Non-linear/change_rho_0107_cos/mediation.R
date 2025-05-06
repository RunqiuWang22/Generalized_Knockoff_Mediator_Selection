args=commandArgs(trailingOnly=TRUE)
iii=as.numeric(args[1])

library(knockoff)
library(glmnet)
####different rho: 
rho_list=seq(0,0.7,by=0.1) #8
rho=rho_list[iii%%length(rho_list)+1]
iii=iii%/%length(rho_list)

###different method
method_list=c("Mymethod","MediationFDR","HIMA","HIMA2","eHIMA","MCP") #6
method=method_list[iii%%length(method_list)+1]
iii=iii%/%length(method_list)

### different a
a_list=seq(0.1,0.5,by=0.1) #5
a=a_list[iii%%length(a_list)+1]
iii=iii%/%length(a_list)

### different b
b_list=seq(0.5,1.5,by=0.2) #6
b=b_list[iii%%length(b_list)+1]
iii=iii%/%length(b_list)

### different p
p_list=c(400) #2
p=p_list[iii%%length(p_list)+1]
iii=iii%/%length(p_list)

###different ss
ss_list=c(0.7) #2
ss=ss_list[iii%%length(ss_list)+1]
iii=iii%/%length(ss_list)

### different standard error
std_list=c(0.4) #1
std=std_list[iii%%length(std_list)+1]
iii=iii%/%length(std_list)

### different setting
setting_list=c("S2") #1
setting=setting_list[iii%%length(setting_list)+1]
iii=iii%/%length(setting_list)

### different sigma
sigma_list=c("CS") #1
sigma=sigma_list[iii%%length(sigma_list)+1]
iii=iii%/%length(sigma_list)

batch=iii+1
rep=20
setwd("./")

source("method.R")
source("SIS.R")
TP = rep(NA,rep)
FDP = rep(NA,rep)
q=0.2

File00=sprintf("TPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s_%d.csv",rho,method,a,b,p,ss,std,setting,sigma,batch)
if (!file.exists(File00)){
  
for (i in 1:rep){
  
  ### generate data
  seeds=1111+i+batch*rep
  set.seed(seeds)
  n=1000
  
  
  ### step 1: data generation
  if (setting=="S2"){
    patha=1:30
    pathb=10:39
  }

  pathi = intersect(patha, pathb)
  if (sigma=="CS"){
    covmat = matrix(rho, nrow = p, ncol = p); diag(covmat) = 1
  }
  
  if (sigma=="AR"){
    covmat = rho^abs(outer(1:p,1:p,"-")); diag(covmat) = 1
  }
  
  data=datagen.nonlinear(n = n, p = p, a = a, b = b, ss=ss, rho = rho,std=std,patha=patha,pathb=pathb,covmat=covmat) 
  
  if (method=="Mymethod1") {
    SIS_select = SIS(X=data$X,Y=data$Y,M=data$M,COV.XM=NULL,Y.family = "gaussian",M.family = "gaussian")
    M1 = data$M[,SIS_select]
    myselect0 = myest(data$X,data$Y,M1,method='Mymethod',q=q)
    myselect = SIS_select[myselect0]
  }
  
  if (method!="Mymethod1") {
    myselect = myest(data$X,data$Y,data$M,V=NULL,method=method,q=q)
  }
  
  File=sprintf("select_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s_%d.csv",rho,method,a,b,p,ss,std,setting,sigma,seeds)
  write.csv(myselect,File)
  
  TP[i]=length(which(myselect%in%c(pathi)))/length(pathi)
  FDP[i]=(length(myselect)-length(which(myselect%in%c(pathi))))/max(1,length(myselect))
  
}


write.csv(TP,File00)

File01=sprintf("FDPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s_%d.csv",rho,method,a,b,p,ss,std,setting,sigma,batch)
write.csv(FDP,File01)
}

