# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 crrSC,
                 cmprsk)
})

# Load data
df <- readRDS("data/readmission90_data.RData")

# Define explanatory variables
explanatory <- c("age_new", "sex", 
                 "country_region", "country_income",
                 "hpd_admreason",
                 "comorbidities_CCI", 
                 "severity_score_scale",
                 "icu_hd_ap", "infection_types",
                 "aci_car", "ent_thir", "ent_car")

# Convert readm_event to factor (required by crr)
df$readm_event <- as.factor(df$readm_event)

# Confirm outcome coding
print(table(df$readm_event))  # ensure "Readmission" is level 1

# Define custom summary function
custom_summary_crrs <- function(model, var_name, var_data) {
  coefs <- model$coef
  var_matrix <- model$var
  se <- sqrt(diag(var_matrix))
  z_values <- coefs / se
  p_values <- 2 * pnorm(-abs(z_values))
  
  # Format p-values
  p_formatted <- ifelse(p_values < 0.001, "<0.001",
                        ifelse(p_values < 0.01, sprintf("%.3f", p_values),
                               sprintf("%.2f", p_values)))
  
  # Exponentiate coefficients
  HR <- exp(coefs)
  CI <- exp(coefs + outer(se, c(-1.96, 1.96)))
  
  summary_table <- data.frame(
    variable = var_name,
    coef = coefs,
    HR = HR,
    lower_95 = CI[,1],
    upper_95 = CI[,2],
    se = se,
    z = z_values,
    p = p_formatted,
    stringsAsFactors = FALSE
  )
  
  # Add factor level names if applicable
  if (is.factor(var_data) || is.character(var_data)) {
    lvls <- levels(as.factor(var_data))
    summary_table$level <- lvls[-1]  # remove reference
  } else {
    summary_table$level <- NA
  }
  
  summary_table <- summary_table %>% select(variable, level, everything())
  return(summary_table)
}

# Run crr() for each variable
results_list <- list()

for (var in explanatory) {
  var_data <- df[[var]]
  
  # Construct design matrix with proper contrast
  covariate_matrix <- model.matrix(~ var_data)[, -1, drop = FALSE]
  
  # Fit Fine-Gray model
  mod <- crr(ftime = df$readm_death_time,
             fstatus = as.numeric(as.character(df$readm_event)),  # convert back to numeric
             cov1 = covariate_matrix,
             failcode = 1,
             cencode = 0)
  
  # Get summary table
  summary_df <- custom_summary_crrs(mod, var, var_data)
  
  results_list[[var]] <- summary_df
}

# Combine results
combined_results <- bind_rows(results_list, .id = NULL)

# -------------------------------------------------------------------
# Plot
hazard_ratio_df <- transform(combined_results, 
                             Hazard_Ratio = sprintf("%.3f", HR), 
                             Lower_CI = sprintf("%.3f", lower_95), 
                             Upper_CI = sprintf("%.3f", upper_95))

hazard_ratio_df$`Hazard ratio (95%CI)` = paste0(hazard_ratio_df$Hazard_Ratio, " [", hazard_ratio_df$Lower_CI, ", ", hazard_ratio_df$Upper_CI, "]", sep = "")

df_used <- hazard_ratio_df %>% 
  select(variable, level, Hazard_Ratio, Lower_CI, Upper_CI, `Hazard ratio (95%CI)`, p)

###
result <- rbind(
  c("Characteristics", NA, NA, NA, "Crude HR (95%CI)", "p.value"),
  c("Age", NA, NA, NA, NA, NA),
  c(NA, df_used[1, c(3:7)]),
  c("Sex", NA, NA, NA, NA, NA),
  c("Female", NA, NA, NA, "Ref", NA),
  c("Male", df_used[2, c(3:7)]),
  c("Region", NA, NA, NA, NA, NA),
  c("Eastern Mediterranean Region", NA, NA, NA, "Ref", NA),
  c("South-East Asian Region", df_used[3, c(3:7)]),
  c("Western Pacific Region", df_used[4, c(3:7)]),
  c("World Bank income status", NA, NA, NA, NA, NA),
  c("High income", NA, NA, NA, "Ref", NA),
  c("Upper middle income", df_used[5, c(3:7)]),
  c("Lower middle income", df_used[6, c(3:7)]),
  c("Primary admission reason", NA, NA, NA, NA, NA),
  c("Infectious disease", NA, NA, NA, "Ref", NA),
  c("Gastrointestinal disorder", df_used[7, c(3:7)]),
  c("Pulmonary disease", df_used[8, c(3:7)]),
  c("Others", df_used[9, c(3:7)]),
  c("Charlson comorbidity index", NA, NA, NA, NA, NA),
  c(NA, df_used[10, c(3:7)]),
  c("Severity score of disease", NA, NA, NA, NA, NA),
  c(NA, df_used[11, c(3:7)]),
  c("Admission to ICU/HD at enrollment", NA, NA, NA, NA, NA),
  c("No", NA, NA, NA, "Ref", NA),
  c("Yes", df_used[12, c(3:7)]),
  c("Infection syndromes", NA, NA, NA, NA, NA),
  c("VAP", NA, NA, NA, "Ref", NA),
  c("Hospital-acquired BSI", df_used[13, c(3:7)]),
  c("Healthcare-associated BSI", df_used[14, c(3:7)]),
  c("Carbapenem-resistant Acinetobacter spp.", NA, NA, NA, NA, NA),
  c("Absent", NA, NA, NA, "Ref", NA),
  c("Present",df_used[15, c(3:7)]),
  c("Third-generation cephalosporin-resistant Enterobacterales", NA, NA, NA, NA, NA),
  c("Absent", NA, NA, NA, "Ref", NA),
  c("Present", df_used[16, c(3:7)]),
  c("Carbapenem-resistant Enterobacterales", NA, NA, NA, NA, NA),
  c("Absent", NA, NA, NA, "Ref", NA),
  c("Present", df_used[17, c(3:7)])
)

#
result <- as.data.frame(result)

result[, 5:6] <- lapply(result[, 5:6], as.character)
result[, 2:4] <- lapply(result[, 2:4], as.numeric)

# Save table
result_save <- result[, c(1, 5, 6)]
saveRDS(result_save, "data/readmission_table_univariable.RData")
###