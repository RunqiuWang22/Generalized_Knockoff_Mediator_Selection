setwd("~/Documents/mediation/data/real_MCI_change_once/")
library(lavaan)
library(ranger)  # For Random Forest
library(dplyr)

for (method in c("Mymethod_RF")) {
  for (Y_column in c("CDRSB","ADAS13", "MMSE")) {
    try({
      ### read the data
      File <- sprintf("final_data_%s.csv", Y_column)
      ddata <- read.csv(File, header=T)
      
      File00 <- sprintf("result/select_%s_%s.csv", method, Y_column)
      select <- read.csv(File00, header=T)
      
      if (method %in% c("HIMA", "HIMA2", "MediationFDR", "eHIMA")) {
        select <- select[,-1]
      }
      
      if (method %in% c("Mymethod_lasso", "Mymethod_RF")) {
        select <- select[,-1]
      }
      
      subdata <- ddata[,names(ddata) %in% c("X", "Y", "age", select)]
      p <- length(select)
      NDE <- NIE <- rep(NA, p)
      final <- matrix(NA, nrow=p, ncol=3)
      
      names(subdata) <- c(paste("M", 1:p, sep=""), "X", "Y", "age")
      
      # Initialize vectors to store results
      alpha_results <- vector("character", p)
      beta_results <- vector("character", p)
      mediation_results <- vector("character", p)
      
      # Function to determine effect direction using partial dependence
      get_effect_direction <- function(model, data, mediator) {
        # Create two datasets with different mediator values
        q1 <- quantile(data[[mediator]], 0.25)
        q3 <- quantile(data[[mediator]], 0.75)
        
        data_low <- data
        data_high <- data
        data_low[[mediator]] <- q1
        data_high[[mediator]] <- q3
        
        pred_low <- predict(model, data = data_low)$predictions
        pred_high <- predict(model, data = data_high)$predictions
        
        # Return direction based on average prediction difference
        return(sign(mean(pred_high - pred_low)))
      }
      
      # Modified function to calculate variable importance and CI
      calculate_rf_importance <- function(data, mediator) {
        tryCatch({
          # Create training data
          train_data <- data.frame(
            M = data[[mediator]],
            X = data$X,
            age = data$age,
            Y = data$Y
          )
          
          # Fit RF model with importance
          rf_model <- ranger(
            Y ~ ., 
            data = train_data,
            num.trees = 1000,
            importance = 'permutation',
            keep.inbag = TRUE
          )
          
          # Get raw importance score and its sign
          raw_imp_score <- importance(rf_model)["M"]
          effect_direction <- get_effect_direction(rf_model, train_data, "M")
          
          # Print diagnostics
          cat("Raw importance score for", mediator, ":", raw_imp_score, "\n")
          cat("Effect direction:", effect_direction, "\n")
          
          # Use absolute value for importance calculations but preserve sign
          imp_score <- abs(raw_imp_score) * effect_direction
          
          # Bootstrap for confidence intervals
          n_boot <- 100
          boot_scores <- numeric(n_boot)
          
          for(i in 1:n_boot) {
            boot_idx <- sample(1:nrow(train_data), replace = TRUE)
            boot_data <- train_data[boot_idx, ]
            
            rf_boot <- ranger(
              Y ~ ., 
              data = boot_data,
              num.trees = 500,
              importance = 'permutation'
            )
            
            boot_score <- importance(rf_boot)["M"]
            boot_direction <- get_effect_direction(rf_boot, boot_data, "M")
            boot_scores[i] <- abs(boot_score) * boot_direction
          }
          
          # Calculate confidence intervals
          ci <- quantile(boot_scores, c(0.025, 0.975))
          
          # Calculate p-value using absolute values
          n_perm <- 100
          perm_scores <- numeric(n_perm)
          
          for(i in 1:n_perm) {
            perm_data <- train_data
            perm_data$M <- sample(perm_data$M)
            
            rf_perm <- ranger(
              Y ~ ., 
              data = perm_data,
              num.trees = 500,
              importance = 'permutation'
            )
            
            perm_scores[i] <- abs(importance(rf_perm)["M"])
          }
          
          # Two-sided p-value
          p_value <- 2 * min(
            mean(perm_scores >= abs(imp_score)),
            mean(perm_scores <= abs(imp_score))
          )
          
          return(list(
            importance = imp_score,
            raw_importance = raw_imp_score,
            direction = effect_direction,
            ci = ci,
            p_value = p_value
          ))
          
        }, error = function(e) {
          warning(paste("Error in RF importance calculation:", e$message))
          return(list(
            importance = NA,
            raw_importance = NA,
            direction = NA,
            ci = c(NA, NA),
            p_value = NA
          ))
        })
      }
      
      # Process each mediator
      for (i in 1:p) {
        m <- paste0("M", i)
        cat("\nProcessing mediator", m, "...\n")
        
        # Path-a: X -> M (linear model)
        tryCatch({
          lm_formula <- as.formula(paste(m, "~ X + age"))
          lm_a <- lm(lm_formula, data = subdata)
          alpha_est <- coef(lm_a)["X"]
          alpha_se <- summary(lm_a)$coefficients["X", "Std. Error"]
          alpha_pvalue <- summary(lm_a)$coefficients["X", "Pr(>|t|)"]
          alpha_results[i] <- paste(round(alpha_est,3), " (", 
                                    round(alpha_est - 1.96 * alpha_se,3), ", ",
                                    round(alpha_est + 1.96 * alpha_se,3), "; ",
                                    round(alpha_pvalue,3), ")", sep="")
        }, error = function(e) {
          alpha_results[i] <- "NA (NA, NA; NA)"
        })
        
        # Path-b: M -> Y (Random Forest importance)
        tryCatch({
          rf_results <- calculate_rf_importance(subdata, m)
          
          if (!is.na(rf_results$importance)) {
            # Include direction in beta results
            beta_results[i] <- paste(round(rf_results$importance, 3), " (",
                                     round(rf_results$ci[1], 3), ", ",
                                     round(rf_results$ci[2], 3), "; ",
                                     round(rf_results$p_value, 3), ")", sep="")
            
            # Calculate mediation effect considering direction
            if (!is.na(alpha_est)) {
              med_effect <- alpha_est * rf_results$importance
              med_se <- sqrt((rf_results$importance^2 * alpha_se^2) + 
                               (alpha_est^2 * ((rf_results$ci[2] - rf_results$ci[1])/(2*1.96))^2))
              mediation_results[i] <- paste(round(med_effect, 3), " (",
                                            round(med_effect - 1.96 * med_se, 3), ", ",
                                            round(med_effect + 1.96 * med_se, 3), "; ",
                                            round(max(alpha_pvalue, rf_results$p_value), 3), ")", sep="")
            } else {
              mediation_results[i] <- "NA (NA, NA; NA)"
            }
          } else {
            beta_results[i] <- "NA (NA, NA; NA)"
            mediation_results[i] <- "NA (NA, NA; NA)"
          }
        }, error = function(e) {
          beta_results[i] <- "NA (NA, NA; NA)"
          mediation_results[i] <- "NA (NA, NA; NA)"
        })
      }
      
      
      # Calculate total and direct effects using RF predictions with all covariates
      rf_total <- ranger(Y ~ X + age, data = subdata)
      rf_direct <- ranger(Y ~ ., data = subdata)
      
      # Function to calculate effects with proper prediction
      calculate_effect <- function(model, data, increment = 1) {
        pred_base <- predict(model, data = data)$predictions
        data_effect <- data
        data_effect$X <- data_effect$X + increment
        pred_effect <- predict(model, data = data_effect)$predictions
        effect <- mean(pred_effect - pred_base)
        effect_se <- sd(pred_effect - pred_base) / sqrt(nrow(data))
        return(list(est = effect, se = effect_se))
      }
      
      # Total effect
      te_results <- calculate_effect(rf_total, subdata)
      TE <- paste(round(te_results$est, 3), " (",
                  round(te_results$est - 1.96 * te_results$se, 3), ", ",
                  round(te_results$est + 1.96 * te_results$se, 3), ")", sep="")
      
      # Direct effect
      de_results <- calculate_effect(rf_direct, subdata)
      NDE <- paste(round(de_results$est, 3), " (",
                   round(de_results$est - 1.96 * de_results$se, 3), ", ",
                   round(de_results$est + 1.96 * de_results$se, 3), ")", sep="")
      
      # Total indirect effect
      nie_est <- te_results$est - de_results$est
      nie_se <- sqrt(te_results$se^2 + de_results$se^2)
      NIE <- paste(round(nie_est, 3), " (",
                   round(nie_est - 1.96 * nie_se, 3), ", ",
                   round(nie_est + 1.96 * nie_se, 3), ")", sep="")
      
      # Create final results dataframe with proper dimensions
      final <- data.frame(
        name = select,
        alpha = alpha_results,
        beta = beta_results,
        Mediation = mediation_results,
        Indirect = rep(NIE, p),
        Direct = rep(NDE, p),
        Total = rep(TE, p)
      )
      
      File0 <- sprintf("result1_%s_%s.csv", method, Y_column)
      write.csv(final, File0)
    })
  }
}