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
                 htmltools,
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

# CRA, 3GCRE, CRE, CRP, VRE, MRSA
df$aci_car <- is_resistant(df$aci_car)
df$ent_thir <- is_resistant(df$ent_thir)
df$ent_car <- is_resistant(df$ent_car)
df$pse_car <- is_resistant(df$pse_car)
df$entc_van <- is_resistant(df$entc_van)
df$sa_meth <- is_resistant(df$sa_meth)

# Delete first28_patient_days <= 0 
df <- df[which(df$first28_patient_days >= 0), ]

# Delete NA
df <- df[complete.cases(df[c("first28_patient_days", 
                             "first28_death")]), ]


# Time and Event
df$first28_death <- ifelse(df$first28_death == 1, "Dead", "Alive")
df$first28_death <- factor(df$first28_death, levels = c("Dead", "Alive"))

df$event <- ifelse(df$first28_death == "Dead", 1, 0)
df$time <- df$first28_patient_days

# Set time = 0 to 1
df$time[which(df$time == 0)] <- 1

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


# Pre data
# Set reference
df$icu_hd_ap <- relevel(factor(df$icu_hd_ap), ref = "No")

# relever for org
columns_to_relevel <- c("aci_car", "ent_thir", "ent_car", 
                        "pse_car", "entc_van", "sa_meth")

df[columns_to_relevel] <- lapply(df[columns_to_relevel], function(x) {
  relevel(factor(x), ref = "Absent")
})
#
saveRDS(df, "data/first28_data.RData")
###
