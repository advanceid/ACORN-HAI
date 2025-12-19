# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 survival,
                 survminer,
                 forestplot)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/first28_data.RData")

#####
# Delete Intercept 
del_row <- c("(Intercept)")

# Model
model <- survreg(formula = Surv(time, event) ~ 
                   age_new + sex + 
                   country_region + country_income +
                   hpd_admreason + comorbidities_CCI + 
                   severity_score_scale + icu_hd_ap +
                   aci_car + ent_thir + ent_car + 
                   infection_types,
                  data = df, dist = "lognormal")

p_model <- data.frame(summary(model)$table)
p_model$ratio <- exp(p_model$Value)
p_model$upper_ci <- exp(p_model$Value + 1.96 * p_model$Std..Error)
p_model$lower_ci <- exp(p_model$Value - 1.96 * p_model$Std..Error)

p_model$p.value <- ifelse(p_model$p < 0.001, "<0.001",
                          ifelse(p_model$p < 0.01, sprintf("%.3f", p_model$p), 
                                 sprintf("%.2f", p_model$p)))
  
  
  ifelse(p_model$p < 0.001, "<0.001", 
                           ifelse(p_model$p < 0.01, "<0.01", 
                                  ifelse(p_model$p < 0.05, "<0.05", sprintf("%.3f", p_model$p))))

p_model <- p_model[-which(rownames(p_model) %in% del_row), ]

ratio_df <- transform(p_model, 
                      Ratio = sprintf("%.3f", ratio), 
                      Lower_CI = sprintf("%.3f", lower_ci), 
                      Upper_CI = sprintf("%.3f", upper_ci))

ratio_df$`ratio (95%CI)` = paste0(ratio_df$Ratio, " [", ratio_df$Lower_CI, ", ", ratio_df$Upper_CI, "]", sep = "")

###
rows <- list(
  c("Characteristics", NA, NA, NA, "Adjusted AF (95%CI)", "P value"),
  c("Age", NA, NA, NA, NA, NA),
  c(NA, ratio_df[1, c(5:7, 12, 8)]),
  c("Sex", NA, NA, NA, NA, NA),
  c("Female", NA, NA, NA, "Ref", NA),
  c("Male", ratio_df[2, c(5:7, 12, 8)]),
  c("Region", NA, NA, NA, NA, NA),
  c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
  c("South-East Asian Region", ratio_df[3, c(5:7, 12, 8)]),
  c("Western Pacific Region", ratio_df[4, c(5:7, 12, 8)]),
  c("World Bank income status", NA, NA, NA, NA, NA),
  c("High income", NA, NA, NA, "Ref", NA),
  c("Upper middle income", ratio_df[5, c(5:7, 12, 8)]),
  c("Lower middle income", ratio_df[6, c(5:7, 12, 8)]),
  c("Primary admission reason", NA, NA, NA, NA, NA),
  c("Infectious disease", NA, NA, NA, "Ref", NA),
  c("Gastrointestinal disorder", ratio_df[7, c(5:7, 12, 8)]),
  c("Pulmonary disease", ratio_df[8, c(5:7, 12, 8)]),
  c("Others", ratio_df[9, c(5:7, 12, 8)]),
  c("Charlson comorbidity index", NA, NA, NA, NA, NA),
  c(NA, ratio_df[10, c(5:7, 12, 8)]),
  c("Severity score of disease", NA, NA, NA, NA, NA),
  c(NA, ratio_df[11, c(5:7, 12, 8)]),
  c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
  c("No", NA, NA, NA, "Ref", NA),
  c("Yes", ratio_df[12, c(5:7, 12, 8)]),
  c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
  c("Absent", NA, NA, NA, "Ref", NA),
  c("Present",ratio_df[13, c(5:7, 12, 8)]),
  c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
  c("Absent", NA, NA, NA, "Ref", NA),
  c("Present", ratio_df[14, c(5:7, 12, 8)]),
  c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
  c("Absent", NA, NA, NA, "Ref", NA),
  c("Present", ratio_df[15, c(5:7, 12, 8)]),
  c("Infection syndromes", NA, NA, NA, NA, NA),
  c("VAP", NA, NA, NA, "Ref", NA),
  c("Hospital-acquired BSI", ratio_df[16, c(5:7, 12, 8)]),
  c("Healthcare-associated BSI", ratio_df[17, c(5:7, 12, 8)])
)


# Bind the rows into a matrix
result <- do.call(rbind, rows)

#
result <- as.data.frame(result)

result[, 5:6] <- lapply(result[, 5:6], as.character)
result[, 2:4] <- lapply(result[, 2:4], as.numeric)
result[is.na(result)] <- NA

# Save table
result_save <- result[, c(1, 5, 6)]
saveRDS(result_save, "data/frist28_table_multi_1_overall.RData")

# -----------------------------------------
# Plot
result_plot <- result

fig <- forestplot(
  result_plot[, c(1, 5, 6)],
  mean = result_plot[, 2],        
  lower = result_plot[, 4],     
  upper = result_plot[, 3],       
  zero = 1,                 
  boxsize = 0.2,              
  graph.pos = 4,             
  hrzl_lines = list(               
    "1" = gpar(lty = 1, lwd = 1.5),  
    "2" = gpar(lty = 2),           
    "40" = gpar(lwd = 1.5, lty = 1, columns = c(1:3)) 
  ),
  graphwidth = unit(.15, "npc"),   
  xlab = "Accelerated factor", 
  xticks = c(seq(0.3, 3, 0.6)), 
  is.summary = c(T, T, F, T, F, F,
                 T, rep(F,3),
                 T, rep(F,3),
                 T, rep(F,4),
                 T, F, T, F,
                 rep(c(T, F, F), 4), 
                 T, rep(F,3)),  
  txt_gp = fpTxtGp(               
    label = gpar(cex = 0.7),      
    ticks = gpar(cex = 0.6),         
    xlab = gpar(cex = 0.6),       
    title = gpar(cex = 0.8)       
  ),
  lwd.zero = 1,                    
  lwd.ci = 1.5,                    
  lwd.xaxis = 1.5,                   
  lty.ci = 1.5,                   
  ci.vertices = T,                 
  ci.vertices.height = 0.1,       
  clip = c(min(result_plot[, 4], na.rm = TRUE), max(result_plot[, 3], na.rm = TRUE)),        
  ineheight = unit(8, 'mm'),      
  line.margin = unit(8, 'mm'),    
  colgap = unit(6, 'mm'),         
  fn.ci_norm = "fpDrawDiamondCI",  
  title = "Multivariable analysis",   
  col = fpColors(               
    box = "blue4",                
    lines = "blue4",               
    zero = "black"                 
  )
)

fig

# Save
pdf(file = "output/pdf/first28_AF.pdf", width = 8, height = 6)
print(fig)
dev.off()
###