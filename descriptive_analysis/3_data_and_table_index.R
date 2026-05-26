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
                 AMR)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df <- readRDS("data/clean_data/baseline_outcomes_index.RData")
ast <- readRDS("data/clean_data/ast_all_index.RData")
hos_reason <- read.xlsx("data/hos_reasons.xlsx")

## Add ast
mdro_positive_levels <- c(
  "Multi-drug-resistant (MDR)",
  "Extensively drug-resistant (XDR)",
  "Pandrug-resistant (PDR)"
)

ast <- ast %>%
  mutate(
    # organism-level: whether this isolate is GNB-MDR
    mdr_gnb = case_when(
      pathogen_group == "GNB" & mdro %in% mdro_positive_levels ~ 1,
      pathogen_group == "GNB" & !is.na(mdro) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(recordid) %>%
  mutate(
    # AMR record-level
    amr_record = case_when(
      any(amr_class == 1, na.rm = TRUE) ~ 1,
      any(amr_class == 0, na.rm = TRUE) ~ 0,
      TRUE ~ NA_real_
    ),
    
    # MDR record-level
    mdr_record = case_when(
      any(mdro %in% mdro_positive_levels, na.rm = TRUE) ~ 1,
      any(!is.na(mdro)) ~ 0,
      TRUE ~ NA_real_
    ),
    
    # MDR-GNB record-level
    mdr_gnb_record = case_when(
      any(mdr_gnb == 1, na.rm = TRUE) ~ 1,
      any(mdr_gnb == 0, na.rm = TRUE) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()


df_add <- ast %>% 
  select(recordid,
         infection_types, 
         pathogen_combined_types,
         pathogen_group_combined,
         org_combined,
         amr = amr_record,
         mdr = mdr_record,
         mdr_gnb = mdr_gnb_record) %>%
  distinct() %>%
  mutate(
    across(c(amr, mdr, mdr_gnb),
           ~ factor(.,
                    levels = c(0, 1),
                    labels = c("-susceptible", "-resistant")))
  ) %>%
  as.data.frame()


###
df_list <- list(df, df_add)

# Convert the `infection_types` column to factor type in all data frames
df_list <- lapply(df_list, function(df) {
  df %>% mutate(infection_types = as.factor(infection_types))
})

df <- reduce(df_list, left_join, by = c("recordid", "infection_types"))

df <- df %>%
  mutate(
    amr = case_when(
      amr == "-resistant" ~ "-resistant",
      amr == "-susceptible" ~ "-susceptible",
      is.na(amr) ~ "-susceptible",
      TRUE ~ as.character(amr)
    ),
    
    mdr = case_when(
      mdr == "-resistant" ~ "-resistant",
      mdr == "-susceptible" ~ "-susceptible",
      is.na(mdr) ~ "-susceptible",
      TRUE ~ as.character(mdr)
    ),
    
    mdr_gnb = case_when(
      mdr_gnb == "-resistant" ~ "-resistant",
      mdr_gnb == "-susceptible" ~ "-susceptible",
      is.na(mdr_gnb) & !is.na(pathogen_group_combined) & grepl("GNB", pathogen_group_combined) ~ "-susceptible",
      TRUE ~ as.character(mdr_gnb)
    )
  ) %>%
  mutate(
    amr = factor(amr, levels = c("-susceptible", "-resistant")),
    mdr = factor(mdr, levels = c("-susceptible", "-resistant")),
    mdr_gnb = factor(mdr_gnb, levels = c("-susceptible", "-resistant"))
  )

###
# Set variables
df_factor <- c("age_group", "country", "country_region", "country_income",
               "infection_types", "sex", "siteid",
               "hpd_admtype", "hpd_admreason", 
               "pathogen_combined_types",
               "readm_ap", "icu_hd_ap", 
               "mdr_gnb", "mdr", "amr",
               "age_group_new", "hai_surg")

df_num <- c("age_raw", "age_new", "age_display_new", "comorbidities_CCI", 
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

###GNB
sum(grepl("GNB", df$pathogen_group_combined), na.rm = TRUE)

#####
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
  
  age_group_new <- factor(age_group_new, 
                          levels = c(1:6), 
                          labels = c("<1 year",
                                     "1–4 years",
                                     "5–14 years",
                                     "15–49 years",
                                     "50–69 years",
                                     "≥70 years")),

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
  hai_surg <- factor(hai_surg, levels = c(1, 0), labels = c("Yes", "No")),
  
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
  age_group_new = "Age group",
  age_display = "Age",
  country_region = "Region",
  country_income = "World Bank income status",
  
  hai_surg = "Surgery before infection onset",
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


# ----------
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
  filter(!is.na(infection_types)) %>%
  mutate(
    infection_types = factor(infection_types, levels = inf_levels),
    age_group_new = factor(age_group_new,
                           levels = c("<1 year",
                                      "1–4 years",
                                      "5–14 years",
                                      "15–49 years",
                                      "50–69 years",
                                      "≥70 years"))
  ) %>%
  group_by(infection_types, age_group_new) %>%
  summarise(
    n = n(),
    med = median(age_display_new, na.rm = TRUE),
    q25 = quantile(age_display_new, 0.25, na.rm = TRUE),
    q75 = quantile(age_display_new, 0.75, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  mutate(
    total = sum(n),
    pct = sprintf("%.1f", 100 * n / total),
    med_s = fmt_age(med),
    q25_s = fmt_age(q25),
    q75_s = fmt_age(q75),
    stat = paste0(
      med_s, " [", q25_s, ", ", q75_s, "]; ",
      n, " (", pct, "%)"
    ),
    stat = gsub("\\b(\\d+)\\.00\\b", "\\1", stat)
  ) %>%
  ungroup() %>%
  select(infection_types, age_group_new, stat) %>%
  pivot_wider(names_from = infection_types, values_from = stat) %>%
  arrange(age_group_new)

# as_tibble()
age_header <- tibble(
  label  = "Age group (years)",
  stat_1 = NA_character_,
  stat_2 = NA_character_,
  stat_3 = NA_character_
)

age_levels <- age_inf_summary %>%
  rename(label = age_group_new) %>%
  mutate(label = paste0("  ", label)) %>%  
  rename(
    stat_1 = `VAP`,
    stat_2 = `Hospital-acquired BSI`,
    stat_3 = `Healthcare-associated BSI`
  ) %>%
  select(label, stat_1, stat_2, stat_3)

age_rows <- bind_rows(age_header, age_levels)

# Overall × age_group
age_overall <- df %>%
  filter(!is.na(age_group_new), !is.na(age_new)) %>%
  mutate(
    age_group_new = factor(age_group_new,
                           levels = c("<1 year",
                                      "1–4 years",
                                      "5–14 years",
                                      "15–49 years",
                                      "50–69 years",
                                      "≥70 years"))
  ) %>%
  group_by(age_group_new) %>%
  summarise(
    n = n(),
    med = median(age_new, na.rm = TRUE),
    q25 = quantile(age_new, 0.25, na.rm = TRUE),
    q75 = quantile(age_new, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total = sum(n),
    pct = sprintf("%.1f", 100 * n / total),
    med_s = fmt_age(med),
    q25_s = fmt_age(q25),
    q75_s = fmt_age(q75),
    Overall = gsub(
      "\\b(\\d+)\\.00\\b", "\\1",
      paste0(
        med_s, " [", q25_s, ", ", q75_s, "]; ",
        n, " (", pct, "%)"
      )
    )
  ) %>%
  select(age_group_new, Overall)


# Other variables
df_export <- df %>%
  select(sex, 
         country_income, hpd_admreason, 
         comorbidities_CCI, 
         severity_score_scale, icu_hd_ap, 
         hai_surg,
         pathogen_combined_types, infection_types) %>%
  tbl_summary(
    by = infection_types,
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2",
      severity_score_scale ~ "continuous2",
      icu_hd_ap ~ "categorical",
      hai_surg ~ "categorical"
    ),
    statistic = list(
      comorbidities_CCI ~ "{median} [{p25}, {p75}]",
      severity_score_scale ~ "{median} [{p25}, {p75}]",
      icu_hd_ap ~ "{n} ({p}%)",
      hai_surg ~ "{n} ({p}%)"
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
         hai_surg,
         pathogen_combined_types, infection_types) %>%
  tbl_summary(
    by = infection_types,
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2",
      severity_score_scale ~ "continuous2",
      icu_hd_ap ~ "categorical",
      hai_surg ~ "categorical"
    ),
    statistic = list(
      comorbidities_CCI ~ "{median} [{p25}, {p75}]",
      severity_score_scale ~ "{median} [{p25}, {p75}]",
      icu_hd_ap ~ "{n} ({p}%)",
      hai_surg ~ "{n} ({p}%)"
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

# ----------------------------------------------------------
# Data for analysis
low_freq_admreason_analysis <- admreason_counts %>%
  filter(percent < 0.05, !is.na(hpd_admreason)) %>%
  pull(hpd_admreason)  

df_analysis <- df %>%
  mutate(hpd_admreason = as.character(hpd_admreason), 
         hpd_admreason = ifelse(hpd_admreason %in% low_freq_admreason_analysis, "Others", hpd_admreason))

# Set hospital reasons levels
df_analysis$hpd_admreason <- factor(df_analysis$hpd_admreason, levels = c("Infectious disease", "Gastrointestinal disorder", "Pulmonary disease", "Others"))

saveRDS(df_analysis, "data/clean_data/data_table_index_new.RData")

# =============================
# Overall table
df_both_delete <- df$recordid[grep("(_BSI)$", df$recordid)] 

df_remain_bsi <- df[which(df$recordid %in% df_both_delete),]

table(df_remain_bsi$infection_types)

df_remain <- df[-which(df$recordid %in% df_both_delete),]
df_remain$age_group <- factor(df_remain$age_group, 
                              levels = c("Adult (Age >= 18 years)",
                                         "Child (Age 1 month - 17 years)",
                                         "Neonate (Age < 28 days)"))

#
df_remain_table <- df_remain %>%
  select(age_new, age_group_new, sex, 
         country_region, country_income, hpd_admreason, 
         comorbidities_CCI, severity_score_scale, hai_surg,
         icu_hd_ap, readm_ap, first_los, first28_icu, first28_mv, 
         pathogen_combined_types, ent_car_3gcr, infection_types) %>%
  tbl_summary(
    missing = "ifany",
    type = list(
      all_continuous() ~ "continuous2",
      severity_score_scale ~ "continuous2",
      icu_hd_ap ~ "categorical",
      readm_ap ~ "categorical",
      hai_surg ~ "categorical"
    ),
    statistic = list(
      age_new ~ "{mean} ± {sd}",
      comorbidities_CCI ~ "{median} [{p25}, {p75}]",
      severity_score_scale ~ "{median} [{p25}, {p75}]",
      first_los ~ "{median} [{p25}, {p75}]",
      first28_icu ~ "{median} [{p25}, {p75}]",
      first28_mv ~ "{median} [{p25}, {p75}]",
      age_group_new ~ "{n} ({p}%)",
      ent_car_3gcr ~ "{n} ({p}%)",
      icu_hd_ap ~ "{n} ({p}%)", 
      readm_ap ~ "{n} ({p}%)",
      hai_surg ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 0,
    missing_text = "Missing"
  ) 

df_remain_table

df_remain_stratified <- df_remain %>%
  mutate(age_group = fct_na_value_to_level(age_group, level = "Missing")) %>%
  select(age_new, age_group_new, sex, 
         country_region, country_income, hpd_admreason, 
         comorbidities_CCI, severity_score_scale, 
         icu_hd_ap, readm_ap, first_los, first28_icu, first28_mv, 
         pathogen_combined_types, ent_car_3gcr) %>%
  tbl_summary(
    by = age_group_new,
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
  filter(!is.na(age_group_new), !is.na(age_new)) %>%
  mutate(age_group = factor(age_group_new,
                            levels = c("<1 year",
                                       "1–4 years",
                                       "5–14 years",
                                       "15–49 years",
                                       "50–69 years",
                                       "≥70 years"))) %>%
  group_by(age_group_new) %>%
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
  select(age_group_new, Overall)

# =============================
# Text for age summary
# Adults >=15 years
# Children 1-14 years
# Infants <1 year (reported in months)
# =============================

fmt_num <- function(x, digits = 1) {
  x <- as.numeric(x)
  if (is.na(x)) return(NA_character_)
  if (abs(x - round(x)) < 1e-8) {
    as.character(round(x))
  } else {
    format(round(x, digits), nsmall = digits, trim = TRUE)
  }
}

# use age_display_new if this is your final harmonised age variable in YEARS
df_age <- df %>%
  filter(!is.na(age_group_new), !is.na(age_display_new)) %>%
  mutate(age_display_new = as.numeric(age_display_new))

# -----------------------------
# Adults: >=15 years
# -----------------------------
adult_sum <- df_age %>%
  filter(age_group_new %in% c("15–49 years", "50–69 years", "≥70 years")) %>%
  summarise(
    med = median(age_display_new, na.rm = TRUE),
    q25 = quantile(age_display_new, 0.25, na.rm = TRUE),
    q75 = quantile(age_display_new, 0.75, na.rm = TRUE)
  )

# -----------------------------
# Children: 1–14 years
# -----------------------------
child_sum <- df_age %>%
  filter(age_group_new %in% c("1–4 years", "5–14 years")) %>%
  summarise(
    med = median(age_display_new, na.rm = TRUE),
    q25 = quantile(age_display_new, 0.25, na.rm = TRUE),
    q75 = quantile(age_display_new, 0.75, na.rm = TRUE)
  )

# -----------------------------
# Infants: <1 year
# convert years -> months
# -----------------------------
infant_sum <- df_age %>%
  filter(age_group_new == "<1 year", age_display_new < 1) %>%
  summarise(
    med = median(age_display_new * 12, na.rm = TRUE),
    q25 = quantile(age_display_new * 12, 0.25, na.rm = TRUE),
    q75 = quantile(age_display_new * 12, 0.75, na.rm = TRUE)
  )


# ———————————————————
# Poly
poly <- df_remain[which(df_remain$pathogen_combined_types == "Polymicrobial"),]

poly_vap <- poly[which(poly$infection_types == "VAP"),]
nrow(poly_vap)/nrow(poly)

poly_with_GNB <- poly_vap[grepl("GNB", poly_vap$pathogen_group_combined), ]
nrow(poly_with_GNB)/nrow(poly_vap)

nrow(poly_vap)
nrow(poly)
nrow(poly_with_GNB)

# ------------------------------------------------------------
# 1. Build organism-level dataset linked to episode info
# ------------------------------------------------------------
df_pathogen <- df %>%
  select(recordid, infection_types, country_region, age_group_new, org_combined) %>%
  distinct()

df_pathogen_long <- df_pathogen %>%
  filter(!is.na(org_combined)) %>%
  mutate(org_combined = str_split(org_combined, pattern = ",\\s*")) %>%
  unnest(org_combined) %>%
  mutate(
    org_combined = str_trim(org_combined),
    org_combined = na_if(org_combined, "")
  ) %>%
  filter(!is.na(org_combined))

df_pathogen_long2 <- df_pathogen_long %>%
  mutate(org_key = str_to_lower(str_squish(org_combined)))

# ------------------------------------------------------------
# 2. Use the already-derived MDRO result to define MDR
#    (do NOT redefine MDR using >=3 resistant antibiotics)
# ------------------------------------------------------------
ast_long2 <- ast %>%
  mutate(
    org_key = str_to_lower(str_squish(org_names_all)),
    mdr_flag = case_when(
      mdro %in% mdro_positive_levels ~ TRUE,
      mdro == "Negative" ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  transmute(
    recordid,
    org_key,
    pathogen_group,
    ris_Carbapenem,
    `ris_Third-generation cephalosporin`,
    Methicillin,
    Vancomycin,
    mdr_flag
  ) %>%
  distinct()

df_pathogen_long_res <- df_pathogen_long2 %>%
  left_join(ast_long2, by = c("recordid", "org_key")) %>%
  select(-org_key) %>%
  mutate(
    mdr_flag = replace_na(mdr_flag, FALSE)
  )

# ------------------------------------------------------------
# 3. Collapse infection types to VAP / BSI
# ------------------------------------------------------------
df_pathogen_long_res <- df_pathogen_long_res %>%
  mutate(
    infection_type2 = ifelse(
      infection_types %in% c("Hospital-acquired BSI", "Healthcare-associated BSI"),
      "BSI",
      as.character(infection_types)
    ),
    infection_type2 = factor(infection_type2, levels = c("VAP", "BSI")),
    age_group_3 = case_when(
      age_group_new == "<1 year" ~ "<1 year",
      age_group_new %in% c("1–4 years", "5–14 years") ~ "1–14 years",
      age_group_new %in% c("15–49 years", "50–69 years", "≥70 years") ~ "≥15 years",
      TRUE ~ NA_character_
    ),
    age_group_3 = factor(age_group_3, levels = c("<1 year", "1–14 years", "≥15 years"))
  )

# ------------------------------------------------------------
# 4. Episode denominator by infection type
# ------------------------------------------------------------
den_episode_inf <- df %>%
  mutate(
    infection_type2 = ifelse(
      infection_types %in% c("Hospital-acquired BSI", "Healthcare-associated BSI"),
      "BSI",
      as.character(infection_types)
    )
  ) %>%
  distinct(recordid, infection_type2) %>%
  count(infection_type2, name = "total_episode")

den_episode_inf

# ------------------------------------------------------------
# 5. Episode denominator by age group
# ------------------------------------------------------------
den_episode_age3 <- df %>%
  mutate(
    age_group_3 = case_when(
      age_group_new == "<1 year" ~ "<1 year",
      age_group_new %in% c("1–4 years", "5–14 years") ~ "1–14 years",
      age_group_new %in% c("15–49 years", "50–69 years", "≥70 years") ~ "≥15 years",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_group_3)) %>%
  distinct(recordid, age_group_3) %>%
  count(age_group_3, name = "total_episode")

den_episode_age3

# ------------------------------------------------------------
# 6. Define phenotype flags at organism row level
# ------------------------------------------------------------
ent_spec_org <- c(
  "E. coli", "K. pneumoniae", "Klebsiella",
  "Enterobacter", "Serratia", "Proteus", "Morganella",
  "Citrobacter", "Providencia", "Raoultella", "Kluyvera",
  "Cronobacter", "Leclercia", "Pluralibacter",
  "Edwardsiella", "Pantoea", "Escherichia", "Salmonella",
  "Kalamiella"
)

df_episode_flags <- df_pathogen_long_res %>%
  mutate(
    CRA = org_combined == "Acinetobacter spp." & ris_Carbapenem == 1,
    CRE = org_combined %in% ent_spec_org & ris_Carbapenem == 1,
    `3GCRE` = org_combined %in% ent_spec_org & `ris_Third-generation cephalosporin` == 1,
    CRP = org_combined == "Pseudomonas spp." & ris_Carbapenem == 1,
    CRKP = org_combined == "K. pneumoniae" & ris_Carbapenem == 1,
    `3GCR_Ecoli` = org_combined == "E. coli" & `ris_Third-generation cephalosporin` == 1,
    `3GCR_Acinetobacter` = org_combined == "Acinetobacter spp." & `ris_Third-generation cephalosporin` == 1,
    VRE = org_combined == "Enterococcus spp." & Vancomycin == 1,
    MRSA = org_combined == "S. aureus" & Methicillin == 1
  ) %>%
  group_by(recordid, infection_type2, age_group_new, age_group_3) %>%
  summarise(
    CRA = any(CRA, na.rm = TRUE),
    CRE = any(CRE, na.rm = TRUE),
    `3GCRE` = any(`3GCRE`, na.rm = TRUE),
    CRP = any(CRP, na.rm = TRUE),
    CRKP = any(CRKP, na.rm = TRUE),
    `3GCR_Ecoli` = any(`3GCR_Ecoli`, na.rm = TRUE),
    `3GCR_Acinetobacter` = any(`3GCR_Acinetobacter`, na.rm = TRUE),
    VRE = any(VRE, na.rm = TRUE),
    MRSA = any(MRSA, na.rm = TRUE),
    MDR = any(mdr_flag, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. Episode-level summary by infection type
# ------------------------------------------------------------
make_episode_table_inf <- function(data, var_name, label) {
  data %>%
    group_by(infection_type2) %>%
    summarise(
      resistant_n = sum(.data[[var_name]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(den_episode_inf, by = "infection_type2") %>%
    mutate(
      phenotype = label,
      percentage = round(100 * resistant_n / total_episode, 1),
      n_over_total = paste0(resistant_n, "/", total_episode)
    ) %>%
    select(infection_type2, phenotype, resistant_n, total_episode, percentage, n_over_total)
}

table_episode_inf <- bind_rows(
  make_episode_table_inf(df_episode_flags, "CRA", "CRA"),
  make_episode_table_inf(df_episode_flags, "CRKP", "carbapenem-resistant K. pneumoniae"),
  make_episode_table_inf(df_episode_flags, "CRP", "CRP"),
  make_episode_table_inf(df_episode_flags, "3GCR_Ecoli", "third-generation cephalosporin-resistant E. coli"),
  make_episode_table_inf(df_episode_flags, "3GCR_Acinetobacter", "third-generation cephalosporin-resistant Acinetobacter spp."),
  make_episode_table_inf(df_episode_flags, "VRE", "VRE"),
  make_episode_table_inf(df_episode_flags, "MRSA", "MRSA"),
  make_episode_table_inf(df_episode_flags, "MDR", "MDR")
) %>%
  arrange(infection_type2, desc(percentage), desc(resistant_n))

table_episode_inf

# Top patterns for VAP / BSI
table_episode_top_inf <- table_episode_inf %>%
  group_by(infection_type2) %>%
  arrange(desc(percentage), desc(resistant_n), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Save data for treatment
saveRDS(df_episode_flags, "data/clean_data/episode_num_for_treatment.RData")

