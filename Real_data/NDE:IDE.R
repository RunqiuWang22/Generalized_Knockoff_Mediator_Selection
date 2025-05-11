setwd("./")
library(lavaan)
for (method in c("MCP", "HIMA", "HIMA2","MediationFDR","eHIMA","Mymethod_lasso")){
  for (Y_column in c("CDRSB","ADAS13", "MMSE")){ 
    try({
      ### read the data
      File=sprintf("final_data_%s.csv",Y_column)
      ddata = read.csv(File,header=T)
      
      File00=sprintf("result/select_%s_%s.csv",method,Y_column)
      select = read.csv(File00,header=T)
      if (method %in% c("HIMA","HIMA2","MediationFDR","eHIMA")){
        select = select[,-1]
      }
      
      if (method %in% c("Mymethod_lasso","Mymethod_RF")){
        select = select[,-1]
      }
      
      subdata = ddata[,names(ddata) %in% c("X","Y","age",select)]
      p = length(select)
      NDE=NIE=rep(NA,p)
      final = matrix(NA,nrow=p,ncol=3)
      
      names(subdata) = c(paste("M",1:p,sep=""),"X", "Y","V")
      
      ### Y is continuous
      #subdata$Y <- as.ordered(subdata$Y)
      
      generate_model <- function(num_mediators) {
        # Initialize model string
        model <- "  # Direct effect from X to Y \n  Y ~ c*X + d*V \n\n" # d*V added to Y
        
        # Add paths from X and V to mediators
        for (i in 1:num_mediators) {
          model <- paste0(model, "  M", i, " ~ a", i, "*X + e", i, "*V \n") # e_i*V added to M_i
        }
        
        # Add paths from mediators and V to outcome Y
        model <- paste0(model, "\n  # Paths from mediators and V to outcome Y \n")
        for (i in 1:num_mediators) {
          model <- paste0(model, "  Y ~ b", i, "*M", i, " \n")
        }
        
        # Add indirect effects for each mediator
        model <- paste0(model, "\n  # Define indirect effects (mediation effects) \n")
        for (i in 1:num_mediators) {
          model <- paste0(model, "  indirect_M", i, " := a", i, " * b", i, " \n")
        }
        
        # Add total indirect effect
        model <- paste0(model, "\n  # Define the total indirect effect \n  total_indirect := ")
        for (i in 1:num_mediators) {
          if (i == 1) {
            model <- paste0(model, "indirect_M", i)
          } else {
            model <- paste0(model, " + indirect_M", i)
          }
        }
        
        # Add total effect
        model <- paste0(model, "\n\n  # Define the total effect \n  total_effect := c + total_indirect")
        
        return(model)
      }
      
      
      
      model <- generate_model(p)
      
      # Fit the model using lavaan
      fit <- sem(model, data = subdata)
      
      # Summarize the results
      params <- summary(fit, standardized = TRUE, fit.measures = TRUE)
      
      #alpha
      alpha <- list()
      for (m in paste0("M", 1:p,sep="")){
        alpha_est <- params$pe$est[params$pe$lhs == m & params$pe$rhs == "X"] 
        alpha_se <- params$pe$se[params$pe$lhs == m & params$pe$rhs == "X"] 
        alpha_pvalue <- params$pe$pvalue[params$pe$lhs == m & params$pe$rhs == "X"] 
        alpha[[m]] <- paste(round(alpha_est,3)," (",  round(alpha_est-1.96* alpha_se,3),", ", round(alpha_est+1.96* alpha_se,3),"; ",round(alpha_pvalue,3),")" ,sep="")
      }
      
      #beta
      beta <- list()
      for (m in paste0("M", 1:p,sep="")){
        beta_est <- params$pe$est[params$pe$lhs == "Y" & params$pe$rhs == m] 
        beta_se <- params$pe$se[params$pe$lhs == "Y" & params$pe$rhs == m]
        beta_pvalue <- params$pe$pvalue[params$pe$lhs == "Y" & params$pe$rhs == m]
        beta[[m]] <- paste(round(beta_est,3)," (",  round(beta_est-1.96* beta_se,3),", ", round(beta_est+1.96* beta_se,3),"; ",round(beta_pvalue,3),")" ,sep="")
      }
      
      #mediators
      MM <- list()
      for (m in paste0("indirect_M", 1:p,sep="")) {
        M_est <- params$pe$est[params$pe$lhs==m]
        M_se <- params$pe$se[params$pe$lhs==m]
        M_pvalue <- params$pe$pvalue[params$pe$lhs==m]
        MM[[m]] <- paste(round(M_est,3)," (",  round(M_est-1.96* M_se,3),", ", round(M_est+1.96* M_se,3),"; ",round(M_pvalue,3),")" ,sep="")
      }
      
      # direct effects
      NDE_est <- params$pe$est[params$pe$lhs == "Y" & params$pe$rhs == "X"]
      NDE_se <- params$pe$se[params$pe$lhs == "Y" & params$pe$rhs == "X"]
      NDE_pvalue <- params$pe$pvalue[params$pe$lhs == "Y" & params$pe$rhs == "X"]
      NDE <- paste(round(NDE_est,3)," (",  round(NDE_est-1.96* NDE_se,3),", ", round(NDE_est+1.96* NDE_se,3),"; ",round(NDE_pvalue,3),")" ,sep="")
      
      # indirect effects 
      NIE_est <- params$pe$est[params$pe$lhs=="total_indirect"]
      NIE_se <- params$pe$se[params$pe$lhs=="total_indirect"]
      NIE_pvalue <- params$pe$pvalue[params$pe$lhs=="total_indirect"]
      NIE <- paste(round(NIE_est,3)," (",  round(NIE_est-1.96* NIE_se,3),", ", round(NIE_est+1.96* NIE_se,3),"; ",round(NIE_pvalue,3),")" ,sep="")
      
      # total effects 
      TE_est <- params$pe$est[params$pe$lhs=="total_effect"]
      TE_se <- params$pe$se[params$pe$lhs=="total_effect"]
      TE_pvalue <- params$pe$pvalue[params$pe$lhs=="total_effect"]
      TE <- paste(round(TE_est,3)," (",  round(TE_est-1.96* TE_se,3),", ", round(TE_est+1.96* TE_se,3),"; ",round(TE_pvalue,3),")" ,sep="")
      
      final = data.frame(name=select,alpha=unlist(alpha),beta=unlist(beta),Mediation=unlist(MM),Indirect=NIE,Direct=NDE,Total=TE)
      File0=sprintf("result_%s_%s.csv",method,Y_column)
      write.csv(final,File0)
    })
  }
}
