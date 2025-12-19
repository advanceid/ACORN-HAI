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
df_all <- readRDS("data/first28_data.RData")
df_used <- split(df_all, df_all$infection_types)

###
# Multi analysis
result_save <- list()

for (i in 1:3){
  df <- df_used[[i]]
  
  # MODEL 
  if (i == 1) {
    # Model 1
    formula <- Surv(time, event) ~ 
      age_new + sex + 
      country_region + country_income +
      hpd_admreason + comorbidities_CCI + 
      severity_score_scale +
      icu_hd_ap +
      aci_car + ent_thir + ent_car
  } else if (i == 2 || i == 3) {
    # Model 2, 3: delete comorbidities_CCI
    formula <- Surv(time, event) ~ 
      age_new + sex + 
      country_region + country_income +
      hpd_admreason + 
      severity_score_scale +
      icu_hd_ap + 
      aci_car + ent_thir + ent_car
  } 
  
  # Fit the model
  model <- survreg(formula, data = df, dist = "lognormal")
  
  p_model <- data.frame(summary(model)$table)
  p_model$ratio <- exp(p_model$Value)
  p_model$upper_ci <- exp(p_model$Value + 1.96 * p_model$Std..Error)
  p_model$lower_ci <- exp(p_model$Value - 1.96 * p_model$Std..Error)
  
  p_model$p.value <- ifelse(p_model$p < 0.001, "<0.001",
                            ifelse(p_model$p < 0.01, sprintf("%.3f", p_model$p), 
                                   sprintf("%.2f", p_model$p)))
  
  del_row <- c("(Intercept)")
  p_model <- p_model[-which(rownames(p_model) %in% del_row), ]
  
  ratio_df <- transform(p_model, 
                        Ratio = sprintf("%.3f", ratio), 
                        Lower_CI = sprintf("%.3f", lower_ci), 
                        Upper_CI = sprintf("%.3f", upper_ci))
  
  ratio_df$`ratio (95%CI)` = paste0(ratio_df$Ratio, " [", ratio_df$Lower_CI, ", ", ratio_df$Upper_CI, "]", sep = "")
  
  ###
  if (i == 1) {
    rows <- list(
      c("Characteristics", NA, NA, NA, "Adjusted AF (95%CI)", "p.value"),
      c("Age", NA, NA, NA, NA, NA),
      c(NA, ratio_df[1, c(9:12, 8)]),
      c("Sex", NA, NA, NA, NA, NA),
      c("Female", NA, NA, NA, "Ref", NA),
      c("Male", ratio_df[2, c(9:12, 8)]),
      c("Region", NA, NA, NA, NA, NA),
      c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
      c("South-East Asian Region", ratio_df[3, c(9:12, 8)]),
      c("Western Pacific Region", ratio_df[4, c(9:12, 8)]),
      c("World Bank income status", NA, NA, NA, NA, NA),
      c("High income", NA, NA, NA, "Ref", NA),
      c("Upper middle income", ratio_df[5, c(9:12, 8)]),
      c("Lower middle income", ratio_df[6, c(9:12, 8)]),
      c("Primary admission reason", NA, NA, NA, NA, NA),
      c("Infectious disease", NA, NA, NA, "Ref", NA),
      c("Gastrointestinal disorder", ratio_df[7, c(9:12, 8)]),
      c("Pulmonary disease", ratio_df[8, c(9:12, 8)]),
      c("Others", ratio_df[9, c(9:12, 8)]),
      c("Charlson comorbidity index", NA, NA, NA, NA, NA),
      c(NA, ratio_df[10, c(9:12, 8)]),
      c("Severity score of disease", NA, NA, NA, NA, NA),
      c(NA, ratio_df[11, c(9:12, 8)]),
      c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
      c("No", NA, NA, NA, "Ref", NA),
      c("Yes", ratio_df[12, c(9:12, 8)]),
      c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c("Present",ratio_df[13, c(9:12, 8)]),
      c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c("Present", ratio_df[14, c(9:12, 8)]),
      c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c("Present", ratio_df[15, c(9:12, 8)])
    )
  } else if (i == 2 || i == 3) {
    rows <- list(
      c("Characteristics", NA, NA, NA, "Adjusted AF (95%CI)", "p.value"),
      c("Age", NA, NA, NA, NA, NA),
      c(NA, ratio_df[1, c(9:12, 8)]),
      c("Sex", NA, NA, NA, NA, NA),
      c("Female", NA, NA, NA, "Ref", NA),
      c("Male", ratio_df[2, c(9:12, 8)]),
      c("Region", NA, NA, NA, NA, NA),
      c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
      c("South-East Asian Region", ratio_df[3, c(9:12, 8)]),
      c("Western Pacific Region", ratio_df[4, c(9:12, 8)]),
      c("World Bank income status", NA, NA, NA, NA, NA),
      c("High income", NA, NA, NA, "Ref", NA),
      c("Upper middle income", ratio_df[5, c(9:12, 8)]),
      c("Lower middle income", ratio_df[6, c(9:12, 8)]),
      c("Primary admission reason", NA, NA, NA, NA, NA),
      c("Infectious disease", NA, NA, NA, "Ref", NA),
      c("Gastrointestinal disorder", ratio_df[7, c(9:12, 8)]),
      c("Pulmonary disease", ratio_df[8, c(9:12, 8)]),
      c("Others", ratio_df[9, c(9:12, 8)]),
      c("Charlson comorbidity index", NA, NA, NA, NA, NA),
      c(NA, NA, NA, NA, NA, NA),
      c("Severity score of disease", NA, NA, NA, NA, NA),
      c(NA, ratio_df[10, c(9:12, 8)]),
      c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
      c("No", NA, NA, NA, "Ref", NA),
      c("Yes", ratio_df[11, c(9:12, 8)]),
      c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c("Present",ratio_df[12, c(9:12, 8)]),
      c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c("Present", ratio_df[13, c(9:12, 8)]),
      c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
      c("Absent", NA, NA, NA, "Ref", NA),
      c("Present", ratio_df[14, c(9:12, 8)])
    )
  } 
  
  # Bind the rows into a matrix
  result <- do.call(rbind, rows)
  
  #
  result <- as.data.frame(result)
  
  result[, 5:6] <- lapply(result[, 5:6], as.character)
  result[, 2:4] <- lapply(result[, 2:4], as.numeric)
  
  result_save[[i]] <- result[, c(1, 5, 6)]
}

# Save table
saveRDS(result_save, "data/first28_table_multi_1.RData")
