library(dplyr)
library(tidyr)
setwd("/Users/runqiuwang/Documents/mediation/simulation/Mediation_continuous_nonlinear/change_rho_0107_cos/change_rho_0107_cos/")
Allfinal=NULL
for (p in c(400)){
  for (a in seq(0.1,0.5,by=0.1)){
    for (b in seq(0.5,1.5,by=0.2)){
      for (ss in c(0.7)){
        for (std in c(0.4)){
          for (setting in c("S2")){
            for (sigma in c("CS")){
              for (rho in seq(0,0.6,by=0.1)){
                Allresult=NULL
                for (method in c("Mymethod","MediationFDR","HIMA","HIMA2","eHIMA","MCP")){
                  for (batch in seq(1:5)){
                    File00=sprintf("TPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s_%d.csv",rho,method,a,b,p,ss,std,setting,sigma,batch)
                    if (file.exists(File00)){
                      TP = read.csv(File00,header=T)
                      TP = TP[,-1]}
                    else{TP=NA}
                    
                    File01=sprintf("FDPP_%.2f_%s_%.2f_%.2f_%d_%.2f_%.2f_%s_%s_%d.csv",rho,method,a,b,p,ss,std,setting,sigma,batch)
                    if (file.exists(File01)){
                      FDP = read.csv(File01,header=T)
                      FDP = FDP[,-1]
                    }
                    else{FDP=NA}
                    result = cbind(method, TP, FDP)
                    Allresult = rbind(Allresult,result)
                  }
                }
                Allresult=as.data.frame(Allresult)
                final <- Allresult %>%
                  mutate(
                    TP = as.numeric(as.character(TP)),  # Convert TP to numeric
                    FDP = as.numeric(as.character(FDP)) # Convert FDP to numeric
                  ) %>%
                  group_by(method) %>%
                  summarize(
                    mean_TP = mean(TP, na.rm = TRUE),  # Calculate mean, ignoring NA values
                    mean_FDR = mean(FDP, na.rm = TRUE)# Calculate mean, ignoring NA values
                  )
                
                final = cbind(p,rho,a,b,final)
                colnames(final) = c("p","rho","a","b","Method","power", "FDR")
                Allfinal=rbind(Allfinal,final)
               
              }
            }
          }
        }
      }
    }
  }
}

df_Allfinal <- Allfinal %>%
  pivot_wider(
    id_cols = c(p, rho, a, b),
    names_from = Method,
    values_from = c(power, FDR)
  )
write.csv(df_Allfinal,"result.csv")
