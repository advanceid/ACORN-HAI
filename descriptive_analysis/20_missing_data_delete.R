# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, 
                 magrittr,
                 tidyr,
                 openxlsx)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/data_table_index_new.RData")

# Count NA by variable
vars_check <- c("age_new", "sex", "country_region", "country_income",
                "hpd_admreason", "icu_hd_ap", "comorbidities_CCI", "severity_score_scale")

# Exactly which variable(s) are missing for each row
df_missing_detail <- df %>%
  mutate(across(all_of(vars_check), is.na, .names = "isna_{col}")) %>%
  filter(if_any(starts_with("isna_"), ~ .x)) %>%
  rowwise() %>%
  mutate(missing_vars = paste(vars_check[which(c_across(starts_with("isna_")))], collapse = ", ")) %>%
  ungroup()

sapply(df[vars_check], function(x) sum(is.na(x)))

df_missing_detail$severity_score_scale

# Delete missing data
df_cleaned <- df %>% drop_na(age_new, sex, 
                             country_region, country_income,
                             hpd_admreason, icu_hd_ap, 
                             comorbidities_CCI, severity_score_scale)

# Save
saveRDS(df_cleaned, "data/clean data/data_table_index_new_delete.RData") 
###
