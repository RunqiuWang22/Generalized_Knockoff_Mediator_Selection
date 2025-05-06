setwd("~/Documents/mediation/simulation/Mediation_continuous/change_rho_confounding/result")

for (p in c(100,400)){
  for (a in c(0.2)){
    for (b in c(0.3,0.5)){
      for (ss in c(0.7)){
        for (std in c(0.4)){
          for (setting in c("S2")){
            for (sigma in c("CS")){
              Allresult=NULL
              for (rho in seq(0,0.9,by=0.1)){
                for (method in c("Mymethod","MCP", "HIMA", "HIMA2","eHIMA","MediationFDR")){
                  
                  File00=sprintf("TPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s.csv",rho,method,a,b,p,ss,std,setting,sigma)
                  if (file.exists(File00)){
                    TP = read.csv(File00,header=T)
                    TP = TP[,-1]}
                  else{TP=NA}
                  
                  File01=sprintf("FDPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s.csv",rho,method,a,b,p,ss,std,setting,sigma)
                  if (file.exists(File01)){
                    FDP = read.csv(File01,header=T)
                    FDP = FDP[,-1]
                  }
                  else{FDP=NA}
                  
                  
                  FDR = mean(FDP,na.rm=T)
                  power = mean(TP,na.rm=T)
                  final = c(rho,method, FDR, power)
                  Allresult = rbind(Allresult,final)
                  colnames(Allresult) = c("rho","Method","FDR", "power")
                }
              }
              File03=sprintf("result_%.2f_%.2f_%d_%.2f_%.2f_%s_%s.csv",a,b,p,ss,std,setting,sigma)
              write.csv(Allresult,File03)
            }
          }
        }
      }
    }
  }
}

