# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr,
                 dplyr,
                 gt,
                 extrafont,
                 gtsummary,
                 openxlsx)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")

# Load fonts
loadfonts() 

# Create a function to check if a variable's value ends with "-resistant"
is_resistant <- function(x) {
  ifelse(grepl("-resistant$", x), 1, 0)
}

# CRA, 3GCRE, CRE, CRP, VRE, CloxRSA
df$aci_car <- is_resistant(df$aci_car)
df$ent_thir <- is_resistant(df$ent_thir)
df$ent_car <- is_resistant(df$ent_car)
df$pse_car <- is_resistant(df$pse_car)
df$entc_van <- is_resistant(df$entc_van)
df$sa_meth <- is_resistant(df$sa_meth)

#
df$hpd_adm_date <- as.Date(df$hpd_adm_date, format = "%Y-%m-%d")
df$readm_ap_1_2 <- as.Date(df$readm_ap_1_2, format = "%Y/%m/%d")
df$readm_ori <- df$hpd_adm_date + 90

#
df <- df %>%
  mutate(readm_90days = case_when(
    readm_ap == "Yes" & !is.na(readm_ap_1_2) & readm_ori >= readm_ap_1_2 ~ 1,
    TRUE ~ 0
  ))

df$readm_observed_end <- as.Date(ifelse(df$readm_90days == 1, 
                                        df$readm_ap_1_2, 
                                        df$readm_ori), origin = "1970-01-01")

df$readm_90days_death <- ifelse(df$mortality_date <= df$readm_observed_end, 2, 0)
#
df$readm_event <- ifelse(!is.na(df$readm_90days_death) & df$readm_90days_death == 2, 
                         df$readm_90days_death, 
                         coalesce(df$readm_90days, df$readm_90days_death))

df$readm_death_end <- as.Date(ifelse(df$readm_event == 2, 
                                     df$mortality_date, 
                                     df$readm_observed_end), origin = "1970-01-01")

df$readm_death_time <- as.numeric(difftime(df$readm_death_end, df$hpd_adm_date, units = "days"))

# Set factors
df <- within(df, {c(
  aci_car <- factor(aci_car, levels = c(1, 0), 
                    labels = c("Present", "Absent")),
  
  ent_car <- factor(ent_car, levels = c(1, 0), 
                    labels = c("Present", "Absent")),
  
  ent_thir <- factor(ent_thir, levels = c(1, 0), 
                     labels = c("Present", "Absent")),
  
  pse_car <- factor(pse_car, levels = c(1, 0), 
                    labels = c("Present", "Absent")),
  
  entc_van <- factor(entc_van, levels = c(1, 0), 
                     labels = c("Present", "Absent")),
  
  sa_meth <- factor(sa_meth, levels = c(1, 0), 
                    labels = c("Present", "Absent"))
  
)})

# Delete readm_death_time <= 0 
# Set time 0 to 1 when event = 2 (died)
df$readm_death_time[which(df$readm_death_time == 0 & df$readm_event == 2)] <- 1

df <- df[which(df$readm_death_time > 0), ]

# Delete NA
df <- df[complete.cases(df[c("readm_death_time", "readm_event")]), ]

#
df <- df %>%
  mutate(readm_death_trans = case_when(
    readm_event == 0 ~ "Alive",
    readm_event == 1 ~ "Readmission",
    readm_event == 2 ~ "Dead",
    TRUE ~ as.character(readm_event)
  ))

df$readm_death_trans <- factor(df$readm_death_trans,
                               levels = c("Readmission", "Dead", "Alive"))

# Revise labels
set_labels <- function(df, labels) {
  for (col in names(labels)) {
    df <- labelled::set_variable_labels(df, !!sym(col) := labels[[col]])
  }
  return(df)
}

labels <- list(
  sex = "Sex",
  age_new = "Age (years)",
  age_group = "Age group",
  country_region = "Region",
  country_income = "World Bank income status",
  
  hpd_admreason = "Primary admission reason",
  
  length_before_onset = "Time duration between hospital admission and infection onset (days)",
  icu_hd_ap = "Admission to ICU/HD at enrollment",
  readm_ap = "Have readmissions",
  
  first_los = "Length of hospital stay (days)",
  first28_icu = "Length of ICU stay (days)",
  first28_mv = "Duration of mechanical ventilation (days)",
  
  pathogen_combined_types = "Type of pathogens",
  
  comorbidities_CCI = "Charlson comorbidity index",
  sofa_score = "SOFA score",
  severity_score_scale = "Severity score of disease",
  fbis_score = "FBIS score",
  pitt_score = "PITT score",
  qpitt_score = "qPITT score",
  eq_5d_3l = "EQ-5D-3L score",
  
  aci_car = "Carbapenem-resistant Acinetobacter spp.", 
  ent_thir = "Third-generation cephalosporin-resistant Enterobacterales", 
  ent_car = "Carbapenem-resistant Enterobacterales", 
  pse_car = "Carbapenem-resistant Pseudomonas spp.", 
  entc_van = "Vancomycin-resistant Enterococcus spp.",
  sa_meth = "Methicillin-resistant Staphylococcus aureus"
)

df <- set_labels(df, labels)

# Set reference
df$icu_hd_ap <- relevel(factor(df$icu_hd_ap), ref = "No")

# relever for org
columns_to_relevel <- c("aci_car", "ent_thir", "ent_car", 
                        "pse_car", "entc_van", "sa_meth")

df[columns_to_relevel] <- lapply(df[columns_to_relevel], function(x) {
  relevel(factor(x), ref = "Absent")
})

# Save
saveRDS(df, "data/readmission90_data.RData")