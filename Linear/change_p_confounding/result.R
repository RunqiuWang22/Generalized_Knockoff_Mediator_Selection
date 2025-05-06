library(ggplot2)
library(tidyr)
library(dplyr)
setwd("~/Documents/mediation/simulation/Mediation_continuous/confounding/change_p_confounding/result/")

generate_and_save_plot <- function(a, b, rho, ss, std, setting, sigma) {   
  # Set working directory   
  
  # Construct the filename based on input parameters   
  file_name <- sprintf("result_%.2f_%.2f_%2f_%.2f_%.2f_%s_%s.csv", a, b, rho, ss, std, setting, sigma)   
  
  # Read the data   
  data <- read.csv(file_name, header = TRUE)   
  data <- data[,-1]   
  data <- data[order(data$p),]   
  names(data)[4]="Power"
  
  # Prepare data for plotting   
  FDR_wider <- pivot_wider(data[,-4], names_from = Method, values_from = FDR)   
  Power_wider <- pivot_wider(data[,-3], names_from = Method, values_from = Power)   
  
  # Convert to long format   
  df_long <- pivot_longer(data, cols = c(FDR, Power), names_to = "Metric", values_to = "Value")   
  
  # Reorder the levels of Metric to ensure the desired facet order   
  df_long$Metric <- factor(df_long$Metric, levels = c("Power", "FDR"))   
  # Custom colors for the legend
  custom_colors <- c(
    "eHIMA" = "#FF9999",      # Light red
    "HIMA2" = "#99CC99",      # Medium green
    "MediationFDR" = "#6699FF", # Medium blue (darker than MCP)
    "HIMA" = "#FFD700",       # Gold
    "MCP" = "#33CCFF",        # Bright cyan blue (more distinct from MediationFDR)
    "Mymethod" = "#FF66B2"    # Vibrant pink (eye-catching)
  )
  
  
  
  custom_linetypes <- c(
    "eHIMA" = "dotdash",
    "HIMA2" = "dashed",
    "MediationFDR" = "dotted",
    "HIMA" = "twodash",
    "MCP" = "longdash",
    "Mymethod" = "solid"  
  )
  
  plot <- ggplot(df_long, aes(x = p, y = Value, color = Method, linetype = Method)) + 
    geom_line(size = 1) + 
    geom_point(size = 2) + 
    # Add consistent horizontal lines
    geom_hline(data = df_long %>% filter(Metric == "FDR"), aes(yintercept = 0.2), linetype = "dashed", color = "grey", size = 1) + 
    geom_hline(data = df_long %>% filter(Metric == "Power"), aes(yintercept = 0.8), linetype = "dashed", color = "grey", size = 1) + 
    geom_hline(data = df_long %>% filter(Metric == "Power"), aes(yintercept = 0), linetype = "solid", color = "black", size = 1) + 
    # Facet layout
    facet_wrap(~ Metric, scales = "free_y", ncol = 1) + 
    scale_x_continuous(limits = c(100, 500), breaks = seq(100, 500, by = 100)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) + 
    scale_color_manual(values = custom_colors) + 
    scale_linetype_manual(values = custom_linetypes) +
    labs(
      title = bquote(a == .(a) ~"," ~ b == .(b) ~"," ~ rho == .(rho)),
      x = "Dimension p",
      y = NULL,
      color = NULL,
      linetype = NULL
    ) + 
    theme_minimal() + 
    theme(
      text = element_text(size = 12),
      legend.position = "none",
      strip.text = element_text(size = 14),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_blank(),     # Remove additional borders
      axis.line = element_line(color = "black") , # Add black x and y axis lines
      plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
    )
  
  # Save the plot to a PDF file   
  output_file <- sprintf("plot_conf_p_%.2f_%.2f_%.2f_%.2f_%.2f_%s_%s.pdf", a, b, rho, ss, std, setting, sigma)   
  ggsave(filename = output_file, plot = plot, width = 3.5, height = 6, dpi = 300)        
  
  return(output_file) 
}



for (rho in c(0.1,0.2)){
  for (a in c(0.2)){
    for (b in c(0.3,0.5)){
      for (ss in c(0.7)){
        for (std in c(0.4)){
          for (setting in c("S2")){
            for (sigma in c("CS")){
output_file <- generate_and_save_plot(a=a, b = b, rho=rho, ss = ss, std=std, setting=setting,sigma=sigma)
      }
    }
  }
}
}
  }
}