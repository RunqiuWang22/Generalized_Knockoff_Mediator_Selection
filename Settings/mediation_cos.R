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


batch=iii+1
rep=20
setwd("./")

source("method.R")
source("SIS.R")
TP = rep(NA,rep)
FDP = rep(NA,rep)
q=0.2



  
for (i in 1:rep){
  
  ### generate data
  seeds=1111+i+batch*rep
  set.seed(seeds)
  n=1000
  
  
  ### step 1: data generation
  
    patha=1:30
    pathb=10:39
    pathi = intersect(patha, pathb)
    covmat = matrix(rho, nrow = p, ncol = p); diag(covmat) = 1
    data=datagen.cos(n = n, p = p, a = a, b = b, ss=ss, rho = rho,std=std,patha=patha,pathb=pathb,covmat=covmat) 
  

    myselect = myest(data$X,data$Y,data$M,V=NULL,method=method,q=q)
  
  
  File=sprintf("select_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s_%d.csv",rho,method,a,b,p,ss,std,setting,sigma,seeds)
  write.csv(myselect,File)
  
  TP[i]=length(which(myselect%in%c(pathi)))/length(pathi)
  FDP[i]=(length(myselect)-length(which(myselect%in%c(pathi))))/max(1,length(myselect))
  
}

File00=sprintf("TPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%d.csv",rho,method,a,b,p,ss,std,batch)
write.csv(TP,File00)

File01=sprintf("FDPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%d.csv",rho,method,a,b,p,ss,std,batch)
write.csv(FDP,File01)


