# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(openxlsx,
                 magrittr,
                 dplyr,
                 tidyr,
                 stringr,
                 purrr,
                 tibble,
                 labelled,
                 gtsummary,
                 gt,
                 htmltools,
                 extrafont,
                 forcats,
                 openxlsx)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df <- readRDS("data/clean data/baseline_outcomes_index.RData")
ast <- readRDS("data/clean data/ast_all_index.RData")
hos_reason <- read.xlsx("data/hos_reasons.xlsx")

# Add ast
df_add <- ast %>% 
  select(recordid, infection_types, 
         pathogen_combined_types, pathogen_group_combined, org_combined) %>%
  distinct()%>%
  as.data.frame()

# Add MDR GNB and MDR
ast_multi <- ast %>% dplyr::select(dplyr::starts_with("ris_"))
row_counts <- rowSums(ast_multi == 1, na.rm = TRUE)

recordid_mdr <- unique(ast$recordid[row_counts >= 3])
recordid_mdr_gnb <- unique(ast$recordid[row_counts >= 3 & ast$pathogen_group == "GNB"])

df$mdr <- ifelse(df$recordid %in% recordid_mdr, "-resistant", "-susceptible")
df$mdr <- factor(df$mdr, levels = c("-susceptible","-resistant"))

df$mdr_gnb <- ifelse(df$recordid %in% recordid_mdr_gnb, "-resistant", "-susceptible")
df$mdr_gnb <- factor(df$mdr_gnb, levels = c("-susceptible","-resistant"))


recordid_dead <- df$recordid[df$first28_death == 1]
recordid_gnb_all <- intersect(recordid_mdr_gnb, df$recordid)

recordid_common <- intersect(recordid_mdr_gnb, recordid_dead)
length(recordid_common)/length(recordid_gnb_all)

# Add AMR
recordid_amr <- unique(ast$recordid[row_counts >= 1])
df$amr <- ifelse(df$recordid %in% recordid_amr, "-resistant", "-susceptible")
df$amr <- factor(df$amr, levels = c("-susceptible","-resistant"))

# Combine
df_list <- list(df, df_add)

# Convert the `infection_types` column to factor type in all data frames
df_list <- lapply(df_list, function(df) {
  df %>% mutate(infection_types = as.factor(infection_types))
})

df <- reduce(df_list, left_join, by = c("recordid", "infection_types"))

###
# Set variables
df_factor <- c("age_group", "country", "country_region", "country_income",
               "infection_types", "sex", "siteid",
               "hpd_admtype", "hpd_admreason", 
               "pathogen_combined_types",
               "readm_ap", "icu_hd_ap", "mdr_gnb")

df_num <- c("age_raw", "age_new", "age_display", "comorbidities_CCI", 
            "first28_icu", "first28_mv", "first_los",
            "readm_ap_1", "icu_hd_ap_1", "mv_ap_1",
            "severity_score_scale", 
            "fbis_score","eq_5d_3l", "pitt_score", "qpitt_score")

df_factor_index <- which(names(df) %in% df_factor)
df_num_index <- which(names(df) %in% df_num)

#
df <- as.data.frame(df)
for (i in df_factor_index){df[,i] = as.factor(df[,i])}
for (i in df_num_index){df[,i] = as.numeric(df[,i])}

###
# Split Primary admission reason
split_condition <- function(condition) {
  # Find the first " "
  pos <- regexpr(" ", condition)
  # 
  code <- substr(condition, 1, pos - 1)
  description <- substr(condition, pos + 1, nchar(condition))
  return(c(code, description))
}

#
split_hos <- t(sapply(hos_reason$Primary.admission.reason, split_condition))

# 
split_hos_reason <- data.frame(abb = split_hos[, 1], 
                               names = split_hos[, 2], 
                               stringsAsFactors = FALSE)

## Set new levels for hpd_admreason
new_levels <- c(setdiff(levels(df$hpd_admreason), "OTH"), "OTH")
df$hpd_admreason <- factor(df$hpd_admreason, levels = new_levels)
level_admreason <- levels(df$hpd_admreason)

# Match new names
match_indices <- match(level_admreason, split_hos_reason$abb)
#
new_level_admreason <- split_hos_reason$names[match_indices]
#
df$hpd_admreason <- factor(df$hpd_admreason, levels = levels(df$hpd_admreason), labels = seq_along(levels(df$hpd_admreason)))

###
df <- within(df, {c(
  infection_types <- factor(infection_types, 
                            levels = c("VAP",
                                       "Hospital-acquired BSI",
                                       "Healthcare-associated BSI")),
  
  sex <- factor(sex, levels = c("F", "M"), labels = c("Female", "Male")),
  
  age_group <- factor(age_group, 
                      levels = c(1:3), 
                      labels = c("Adult (Age >= 18 years)",
                                 "Child (Age 1 month - 17 years)",
                                 "Neonate (Age < 28 days)")),
  
  country_income <- factor(country_income, 
                           levels = c("High income", "Upper middle income",
                                      "Lower middle income")),
  
  pathogen_combined_types <- factor(pathogen_combined_types,
                                    levels = c("Gram-negative bacteria",
                                               "Gram-positive bacteria",
                                               "Fungi",
                                               "Polymicrobial")),
  
  hpd_admreason <- factor(hpd_admreason,
                          levels = c(1:17),
                          labels = new_level_admreason),
  
  icu_hd_ap <- factor(icu_hd_ap, levels = c(1, 2), labels = c("Yes", "No")),
  readm_ap <- factor(readm_ap, levels = c(1, 2), labels = c("Yes", "No")),
  
  aci_car <- factor(aci_car, levels = c(1, 0), 
                    labels = c("Carbapenem-resistant", "Carbapenem-susceptible")),
  
  ent_car <- factor(ent_car, levels = c(1, 0), 
                    labels = c("Carbapenem-resistant", 
                               "Carbapenem-susceptible")),
  
  ent_thir <- factor(ent_thir, levels = c(1, 0), 
                     labels = c("Third-generation cephalosporin-resistant", 
                                "Third-generation cephalosporin-susceptible")),
  
  ent_car_3gcr <- factor(ent_car_3gcr, levels = c(2, 1, 0), 
                         labels = c("Carbapenem-resistant", 
                                    "Third-generation cephalosporin-resistant", 
                                    "Third-generation cephalosporin-susceptible")),
  
  pse_car <- factor(pse_car, levels = c(1, 0), 
                    labels = c("Carbapenem-resistant", 
                               "Carbapenem-susceptible")),
  
  entc_van <- factor(entc_van, levels = c(1, 0), 
                     labels = c("Vancomycin-resistant", 
                                "Vancomycin-susceptible")),
  
  sa_meth <- factor(sa_meth, levels = c(1, 0), 
                    labels = c("Methicillin-resistant", 
                               "Methicillin-susceptible"))
  
)})

# Combine hospital reasons less than 5% to "Others"
admreason_counts <- df %>%
  dplyr::count(infection_types, hpd_admreason) %>%
  group_by(infection_types) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

low_freq_admreason <- admreason_counts %>%
  group_by(hpd_admreason) %>%
  summarise(all_below_5 = all(percent < 0.05)) %>%
  filter(all_below_5) %>%
  pull(hpd_admreason)  

df <- df %>%
  mutate(hpd_admreason = as.character(hpd_admreason), 
         hpd_admreason = ifelse(hpd_admreason %in% low_freq_admreason, "Others", hpd_admreason))

# Set hospital reasons levels
current_levels <- unique(df$hpd_admreason)

middle_levels <- sort(setdiff(current_levels, c("Infectious disease", "Others")))

df$hpd_admreason <- factor(df$hpd_admreason, 
                           levels = c("Infectious disease", middle_levels, "Others"))

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
  age_display = "Age",
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
  
  aci_car = "Acinetobacter spp.",
  ent_car_3gcr = "Enterobacterales",
  pse_car = "Pseudomonas spp.",
  entc_van = "Enterococcus spp.",
  sa_meth = "Staphylococcus aureus"
)

df <- set_labels(df, labels)

# Filter those were included in both VAP and BSI groups
df_both <- df[grep("(_VAP|_BSI)$", df$recordid), ] %>%
  group_by(infection_types) %>%
  dplyr::summarise(n = n())

df_both

# =============================
# Age group (years)
inf_levels <- c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")

# <1 2f
fmt_age <- function(x, eps = 1e-9) {
  y <- suppressWarnings(as.numeric(x))
  out <- ifelse(
    is.na(y), NA_character_,
    ifelse(
      y < (1 - eps),
      sprintf("%.2f", y),
      as.character(round(y))
    )
  )
  out
}

# infection_types × age_group
age_inf_summary <- df %>%
  filter(!is.na(infection_types), !is.na(age_group), !is.na(age_display)) %>%
  mutate(infection_types = factor(infection_types, levels = inf_levels)) %>%
  group_by(infection_types, age_group) %>%
  summarise(
    med = median(age_display, na.rm = TRUE),
    q25 = quantile(age_display, 0.25, na.rm = TRUE),
    q75 = quantile(age_display, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    med_s  = fmt_age(med),
    q25_s  = fmt_age(q25),
    q75_s  = fmt_age(q75),
    stat   = paste0(med_s, " [", q25_s, ", ", q75_s, "]"),
    stat   = gsub("\\b(\\d+)\\.00\\b", "\\1", stat)
  ) %>%
  select(infection_types, age_group, stat) %>%
  pivot_wider(names_from = infection_types, values_from = stat) %>%
  arrange(factor(age_group, levels = c("Adult (Age >= 18 years)",
                                       "Child (Age 1 month - 17 years)",
                                       "Neonate (Age < 28 days)")))

# as_tibble()
age_header <- tibble(
  label  = "Age group (years)",
  stat_1 = NA_character_,
  stat_2 = NA_character_,
  stat_3 = NA_character_
)

age_levels <- age_inf_summary %>%
  rename(label = age_group) %>%
  mutate(label = paste0("  ", label)) %>%  
  rename(
    stat_1 = `VAP`,
    stat_2 = `Hospital-acquired BSI`,
    stat_3 = `Healthcare-associated BSI`
  ) %>%
  select(label, stat_1, stat_2, stat_3)

age_rows <- bind_rows(age_header, age_levels)

# Overall × age_group
age_overall_summary <- df %>%
  filter(!is.na(age_group), !is.na(age_display)) %>%
  group_by(age_group) %>%
  summarise(
    med = median(age_display, na.rm = TRUE),
    q25 = quantile(age_display, 0.25, na.rm = TRUE),
    q75 = quantile(age_display, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(infection_types = "Overall")


# Other variables
df_export <- df %>%
  select(sex, 
         country_income, hpd_admreason, 
         comorbidities_CCI, 
         severity_score_scale, icu_hd_ap, 
         pathogen_combined_types, infection_types) %>%
  tbl_summary(
    by = infection_types,
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2",
      severity_score_scale ~ "continuous2",
      icu_hd_ap ~ "categorical"
    ),
    statistic = list(
      comorbidities_CCI ~ "{median} [{p25}, {p75}]",
      severity_score_scale ~ "{median} [{p25}, {p75}]",
      icu_hd_ap ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 0,
      all_categorical() ~ c(0, 1) 
    ),
    missing_text = "Missing"
  ) %>%
  as_tibble() %>%
  mutate(across(where(is.character), ~ gsub(",", "", .x))) 

df_export_overall <- df %>%
  select(sex, 
         country_income, hpd_admreason, 
         comorbidities_CCI, 
         severity_score_scale, icu_hd_ap, 
         pathogen_combined_types, infection_types) %>%
  tbl_summary(
    by = infection_types,
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2",
      severity_score_scale ~ "continuous2",
      icu_hd_ap ~ "categorical"
    ),
    statistic = list(
      comorbidities_CCI ~ "{median} [{p25}, {p75}]",
      severity_score_scale ~ "{median} [{p25}, {p75}]",
      icu_hd_ap ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 0,
    missing_text = "Missing"
  ) %>%
  add_overall() %>%
  as_tibble() %>%
  mutate(across(where(is.character), ~ gsub(",", "", .x)))

# Footnotes
notes <- data.frame(
  Notes = c(
    "Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, IQR = Interquartile Range, ICU = Intensive Care Unit, HD = High Dependency.",
    "146 patients had both index VAP and BSI episodes (135 with hospital-acquired BSI and 11 with healthcare-associated BSI).",
    "Standardized severity score included qSOFA for adults, sepsis six recognition features for children, and general WHO severity signs for neonates.",
    "Median [IQR]"
  ),
  stringsAsFactors = FALSE
)

# Save
wb <- createWorkbook()
addWorksheet(wb, "AgeGroup_median_IQR")
addWorksheet(wb, "OtherFeatures")
addWorksheet(wb, "Notes")

writeData(wb, "AgeGroup_median_IQR", age_rows, keepNA = F)
writeData(wb, "OtherFeatures", df_export, keepNA = F)
writeData(wb, "Notes", notes, keepNA = FALSE)
saveWorkbook(wb, "output/table/summary_baseline.xlsx", overwrite = TRUE)

# =============================
# Data for analysis
low_freq_admreason_analysis <- admreason_counts %>%
  filter(percent < 0.05, !is.na(hpd_admreason)) %>%
  pull(hpd_admreason)  

df_analysis <- df %>%
  mutate(hpd_admreason = as.character(hpd_admreason), 
         hpd_admreason = ifelse(hpd_admreason %in% low_freq_admreason_analysis, "Others", hpd_admreason))

# Set hospital reasons levels
df_analysis$hpd_admreason <- factor(df_analysis$hpd_admreason, levels = c("Infectious disease", "Gastrointestinal disorder", "Pulmonary disease", "Others"))

saveRDS(df_analysis, "data/clean data/data_table_index_new.RData")

# =============================
# Overall table
df_both_delete <- df$recordid[grep("(_BSI)$", df$recordid)] 

df_remain_bsi <- df[which(df$recordid %in% df_both_delete),]

table(df_remain_bsi$infection_types)

df_remain <- df[-which(df$recordid %in% df_both_delete),]
df_remain$age_group <- factor(df_remain$age_group, 
                              levels = c("Adult (Age >= 18 years)",
                                         "Child (Age 1 month - 17 years)",
                                         "Neonate (Age <= 28 days)"))

#
df_remain_table <- df_remain %>%
  select(age_new, age_group, sex, 
         country_region, country_income, hpd_admreason, 
         comorbidities_CCI, severity_score_scale, 
         icu_hd_ap, readm_ap, first_los, first28_icu, first28_mv, 
         pathogen_combined_types, ent_car_3gcr, infection_types) %>%
  tbl_summary(
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2",
      severity_score_scale ~ "continuous2",
      icu_hd_ap ~ "categorical",
      readm_ap ~ "categorical"
    ),
    statistic = list(
      age_new ~ "{mean} ± {sd}",
      comorbidities_CCI ~ "{median} [{p25}, {p75}]",
      severity_score_scale ~ "{median} [{p25}, {p75}]",
      first_los ~ "{median} [{p25}, {p75}]",
      first28_icu ~ "{median} [{p25}, {p75}]",
      first28_mv ~ "{median} [{p25}, {p75}]",
      age_group ~ "{n} ({p}%)",
      ent_car_3gcr ~ "{n} ({p}%)",
      icu_hd_ap ~ "{n} ({p}%)", 
      readm_ap ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 0,
    missing_text = "Missing"
  ) 

df_remain_table

df_remain_stratified <- df_remain %>%
  mutate(age_group = fct_na_value_to_level(age_group, level = "Missing")) %>%
  select(age_new, age_group, sex, 
         country_region, country_income, hpd_admreason, 
         comorbidities_CCI, severity_score_scale, 
         icu_hd_ap, readm_ap, first_los, first28_icu, first28_mv, 
         pathogen_combined_types, ent_car_3gcr) %>%
  tbl_summary(
    by = age_group,
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2"
    ),
    statistic = list(
      all_continuous() ~ "{median} [{p25}, {p75}]",
      icu_hd_ap ~ "{n} ({p}%)", 
      readm_ap ~ "{n} ({p}%)",
      ent_car_3gcr ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 0,
    missing_text = "Missing"
  ) %>%
  add_overall()


df_remain_stratified


age_overall <- df %>%
  filter(!is.na(age_group), !is.na(age_new)) %>%
  mutate(age_group = factor(age_group,
                            levels = c("Adult (Age >= 18 years)",
                                       "Child (Age 1 month - 17 years)",
                                       "Neonate (Age <= 28 days)"))) %>%
  group_by(age_group) %>%
  summarise(
    med = median(age_new, na.rm = TRUE),
    q25 = quantile(age_new, 0.25, na.rm = TRUE),
    q75 = quantile(age_new, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    med_s = fmt_age(med),
    q25_s = fmt_age(q25),
    q75_s = fmt_age(q75),
    Overall = gsub("\\b(\\d+)\\.00\\b", "\\1",
                   paste0(med_s, " [", q25_s, ", ", q75_s, "]"))
  ) %>%
  select(age_group, Overall)

# =============================
# Poly
poly <- df_remain[which(df_remain$pathogen_combined_types == "Polymicrobial"),]

poly_vap <- poly[which(poly$infection_types == "VAP"),]
nrow(poly_vap)/nrow(poly)

poly_with_GNB <- poly_vap[grepl("GNB", poly_vap$pathogen_group_combined), ]
nrow(poly_with_GNB)/nrow(poly_vap)

nrow(poly_vap)
nrow(poly)
nrow(poly_with_GNB)
