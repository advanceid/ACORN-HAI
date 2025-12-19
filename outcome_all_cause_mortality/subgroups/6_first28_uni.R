# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 ggplot2,
                 MASS,
                 car,
                 gt,
                 purrr,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df_all <- readRDS("data/first28_data.RData")
df_used <- split(df_all, df_all$infection_types)

# Univariable analysis
result_save <- list()

for (j in 1:3) {
  df <- df_used[[j]]
  
  # potential variables, delete income
  variables <- c("age_new", "sex", 
                 "country_region", "country_income",
                 "hpd_admreason", "comorbidities_CCI", 
                 "severity_score_scale",
                 "icu_hd_ap", 
                 "aci_car", "ent_thir", "ent_car")
  
  # Initialize the data frame to store results
  p_values <- data.frame(p.value = character(), ratio = numeric(), 
                         lower_ci = numeric(), upper_ci = numeric(), 
                         stringsAsFactors = FALSE)
  
  for (var in variables) {
    formula <- as.formula(paste("Surv(time, event) ~", var))
    univariate_model <- survreg(formula, data = df, dist = "lognormal")
    model_summary <- summary(univariate_model)
    
    # Skip Intercept and scale
    variable_names <- rownames(model_summary$table)
    values <- model_summary$table[, 1]
    se_values <- model_summary$table[, 2]
    p_values_for_var <- model_summary$table[, 4]
    
    for (i in seq_along(variable_names)) {
      if (!grepl("Intercept|Log\\(scale\\)", variable_names[i])) {
        p_value <- p_values_for_var[i]
        value <- values[i]
        se <- se_values[i]
        
        ratio <- exp(value)
        upper_ci <- exp(value + 1.96 * se)
        lower_ci <- exp(value - 1.96 * se)
        
        significance <- ifelse(p_value < 0.001, "<0.001",
                               ifelse(p_value < 0.01, sprintf("%.3f", p_value), 
                                      sprintf("%.2f", p_value)))
        
        
        p_values <- rbind(p_values, data.frame(p.value = significance, ratio = ratio, lower_ci = lower_ci, upper_ci = upper_ci))
      }
    }
  }
  
  # Display the final data frame
  p_values
  
  # -------------------------------------------------------------------
  ratio_df <- transform(p_values, 
                        Ratio = sprintf("%.3f", ratio), 
                        Lower_CI = sprintf("%.3f", lower_ci), 
                        Upper_CI = sprintf("%.3f", upper_ci))
  
  ratio_df$`ratio (95%CI)` = paste0(ratio_df$Ratio, " [", ratio_df$Lower_CI, ", ", ratio_df$Upper_CI, "]", sep = "")
  
  ###
  rows <- list(
    c("Characteristics", NA, NA, NA, "Crude AF (95%CI)", "p.value"),
    c("Age", NA, NA, NA, NA, NA),
    c(NA, ratio_df[1, c(5:8, 1)]),
    c("Sex", NA, NA, NA, NA, NA),
    c("Female", NA, NA, NA, "Ref", NA),
    c("Male", ratio_df[2, c(5:8, 1)]),
    c("Region", NA, NA, NA, NA, NA),
    c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
    c("South-East Asian Region", ratio_df[3, c(5:8, 1)]),
    c("Western Pacific Region", ratio_df[4, c(5:8, 1)]),
    c("World Bank income status", NA, NA, NA, NA, NA),
    c("High income", NA, NA, NA, "Ref", NA),
    c("Upper middle income", ratio_df[5, c(5:8, 1)]),
    c("Lower middle income", ratio_df[6, c(5:8, 1)]),
    c("Primary admission reason", NA, NA, NA, NA, NA),
    c("Infectious disease", NA, NA, NA, "Ref", NA),
    c("Gastrointestinal disorder", ratio_df[7, c(5:8, 1)]),
    c("Pulmonary disease", ratio_df[8, c(5:8, 1)]),
    c("Others", ratio_df[9, c(5:8, 1)]),
    c("Charlson comorbidity index", NA, NA, NA, NA, NA),
    c(NA, ratio_df[10, c(5:8, 1)]),
    c("Severity score of disease", NA, NA, NA, NA, NA),
    c(NA, ratio_df[11, c(5:8, 1)]),
    c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
    c("No", NA, NA, NA, "Ref", NA),
    c("Yes", ratio_df[12, c(5:8, 1)]),
    c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
    c("Absent", NA, NA, NA, "Ref", NA),
    c("Present",ratio_df[13, c(5:8, 1)]),
    c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
    c("Absent", NA, NA, NA, "Ref", NA),
    c("Present", ratio_df[14, c(5:8, 1)]),
    c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
    c("Absent", NA, NA, NA, "Ref", NA),
    c("Present", ratio_df[15, c(5:8, 1)])
  )
  
  # Bind the rows into a matrix
  result <- do.call(rbind, rows)
  result <- as.data.frame(result)
  
  result[, 5:6] <- lapply(result[, 5:6], as.character)
  result[, 2:4] <- lapply(result[, 2:4], as.numeric)
  
  # Save table
  result_save[[j]] <- result[, c(1, 5, 6)]
}

# Save
saveRDS(result_save, "data/first28_table_univariable.RData")
