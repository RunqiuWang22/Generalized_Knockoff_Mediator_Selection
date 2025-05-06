library(dplyr)
setwd("~/Documents/mediation/simulation/Mediation_continuous_nonlinear/change_rho_0107_cos/")
Allfinal=NULL
for (p in c(400)){
 for (a in c(0.2,0.5,1)){
  for (b in c(0.5,0.7,1)){
    for (ss in c(0.7)){
     for (std in c(0.4)){
       for (setting in c("S2")){
         for (sigma in c("CS")){
          for (rho in c(0,0.3,0.5,0.7)){
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
    mean_FDR = mean(FDP, na.rm = TRUE) # Calculate mean, ignoring NA values
  )

 final = cbind(p,a,b,rho,final)
 colnames(final) = c("p","a","b","rho","Method","power", "FDR")
 Allfinal=rbind(Allfinal,final)
 # File03=sprintf("result_%.2f_%.2f_%d_%.2f_%.2f_%s_%s.csv",a,b,p,ss,std,setting,sigma)
 # write.csv(final,File03)
     }
    }
   }
  }
}
}
}
}

### transfer
Allfinal_power = subset(Allfinal,select=c("p","rho","a" ,"b","Method","power"))
power <-  Allfinal_power %>%
  pivot_wider(
    names_from = Method,     # The column whose values become new column names
    values_from = power      # The column whose values fill the new columns
  )
power=power[order(power$rho),]
power=cbind("power",power)
names(power)[1]="result"

Allfinal_FDR = subset(Allfinal,select=c("p","rho","a" ,"b","Method","FDR"))
FDR<-  Allfinal_FDR %>%
  pivot_wider(
    names_from = Method,     # The column whose values become new column names
    values_from = FDR      # The column whose values fill the new columns
  )
FDR=FDR[order(FDR$rho),]
FDR=cbind("FDR",FDR)
names(FDR)[1]="result"
result=rbind(FDR,power)

desired_order <- c("result","p","rho","a","b" ,"HIMA", "HIMA2", "eHIMA", "MCP", "MediationFDR", "Mymethod")

# Reorder columns
result_reordered <- result %>%
  select(all_of(desired_order))


write.csv(result_reordered,"result.csv")


