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

# CRA, 3GCRE, CRE, CRP, VRE, MRSA
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

# ------------------
# Check data
df_check <- df[which(df$readm_death_time <= 0), ] %>%
  select(recordid, hpd_adm_date, readm_ap_1_2, mortality_date, readm_death_time)

write.xlsx(df_check, "readmission_issues.xlsx")


# ------------------

# Delete readm_death_time <= 0 
# Set time 0 to 1 when event = 2 (died)
check_row1 <- which(df$readm_death_time <= 0)
check_row2 <- which(df$readm_death_time == 0 & df$readm_event == 2)
check_row <- setdiff(check_row1, check_row2)

df_check_new <- df[check_row, ]

#
df$readm_death_time[which(df$readm_death_time == 0 & df$readm_event == 2)] <- 1

df_dele <- subset(df, readm_death_time <= 0 | is.na(readm_death_time))


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

# Create the summary table
table <- df %>%
  select(age_new, sex, 
         country_region, country_income, hpd_admreason, 
         comorbidities_CCI, severity_score_scale, 
         icu_hd_ap, aci_car, ent_thir, ent_car, 
         pse_car, entc_van, sa_meth,
         infection_types, readm_death_trans) %>%
  tbl_strata(
    strata = readm_death_trans,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = infection_types,
                  missing = "ifany",
                  type = list(
                    all_continuous() ~ "continuous2",
                    severity_score_scale ~ "continuous2",
                    comorbidities_CCI ~ "continuous2",
                    age_new ~ "continuous2",
                    icu_hd_ap ~ "categorical",
                    aci_car ~ "categorical", 
                    ent_thir ~ "categorical",
                    ent_car ~ "categorical", 
                    pse_car ~ "categorical", 
                    entc_van ~ "categorical",
                    sa_meth ~ "categorical"
                  ),
                  statistic = list(
                    age_new ~ "{median} [{p25}, {p75}]",
                    comorbidities_CCI ~ "{median} [{p25}, {p75}]",
                    severity_score_scale ~ "{median} [{p25}, {p75}]"
                  ),
                  digits = all_continuous() ~ 0,
                  missing_text = "Missing") %>%
      add_overall(),
    .header = "**{strata}, N = {format(n, big.mark = '')}**"
  ) %>%
  modify_fmt_fun(
    update = all_stat_cols() ~ function(x) {
      gsub("(?<=\\d)[,，](?=\\d{3}(\\D|$))", "", x, perl = TRUE)
    }
  ) %>%
  modify_header(
    update = list(
      label ~ "**Characteristics**",
      all_stat_cols() ~ gt::html("**{level}**<br>N = {format(n, big.mark = '')}")
    )
  ) %>%
  bold_labels() %>%
  modify_footnote(update = everything() ~ NA) %>%
  modify_table_styling(
    columns = label,
    rows = label == "Severity score of disease",
    footnote = "Standardized severity score included qSOFA for adults, sepsis six recognition features for children, and general WHO severity signs for neonates.") %>% 
  as_gt() %>%
  tab_options(footnotes.marks = c("†")) %>%
  tab_source_note(
    source_note = "Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, ICU = Intensive Care Unit, HD = High Dependency, IQR = Interquartile Range (25th percentile [Q1] to 75th percentile [Q3])."
  ) %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1:5)
  ) 


# Customize the headers to split lines
gt_table <- table %>%
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2),
    table_body.hlines.color = "white",
    table_body.hlines.width = px(0),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(0),
    data_row.padding = px(0),
    footnotes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(
      rows = nrow(table$`_data`)
    )
  ) %>%
  cols_width(
    1:5 ~ px(195),
    c(6, 7, 11, 12, 16, 17) ~ px(70),
    c(8, 13, 18) ~ px(110),
    c(9, 14, 19) ~ px(130)
    
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10)
    ),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10)
    ),
    locations = cells_footnotes()
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10),
      weight = "bold"
    ),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Times New Roman"),
      size = px(10)
    ),
    locations = cells_source_notes()
  ) %>%
  text_transform(
    locations = cells_body(
      columns = label,
      rows = label == "Carbapenem-resistant Acinetobacter spp."
    ),
    fn = function(x) html("<b>Carbapenem-resistant <i>Acinetobacter</i> spp.</b>")
  ) %>%
  text_transform(
    locations = cells_body(
      columns = label,
      rows = label == "Carbapenem-resistant Pseudomonas spp."
    ),
    fn = function(x) html("<b>Carbapenem-resistant <i>Pseudomonas</i> spp.</b>")
  ) %>%
  text_transform(
    locations = cells_body(
      columns = label,
      rows = label == "Vancomycin-resistant Enterococcus spp."
    ),
    fn = function(x) html("<b>Vancomycin-resistant <i>Enterococcus</i> spp.</b>")
  ) %>%
  text_transform(
    locations = cells_body(
      columns = label,
      rows = label == "Methicillin-resistant Staphylococcus aureus"
    ),
    fn = function(x) html("<b>Methicillin-resistant <i>Staphylococcus aureus</i></b>")
  )

# Print the table
gt_table

# Save
gtsave(gt_table, filename = "output/table/table_summary.html")

# ----------------------------------------------------------------------------
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


# ===========================
# Patient-level 90-day readmission and mortality crude rates
# ===========================

library(dplyr)
library(stringr)
library(scales)

# 1) Clean IDs: drop infection-episode suffixes
df_clean_id <- df %>%
  mutate(recordid_base = str_remove(recordid, "_(VAP|BSI)$"))

# 2) Collapse to patient level: flags for readmission and death
df_patient_flags <- df_clean_id %>%
  group_by(recordid_base) %>%
  summarise(
    ever_readm = as.integer(any(readm_event == 1, na.rm = TRUE)),  # readmission
    ever_death = as.integer(any(readm_event == 2, na.rm = TRUE)),  # death
    .groups = "drop"
  ) %>%
  mutate(
    status_patient = case_when(
      ever_death == 1 ~ "Dead",
      ever_readm == 1 ~ "Readmission",
      TRUE            ~ "Alive"
    )
  )

# 3) Crude readmission rate
crude_readmission_unique <- df_patient_flags %>%
  summarise(
    patients   = n(),
    readmitted = sum(ever_readm, na.rm = TRUE),
    crude      = readmitted / patients,
    crude_pct  = percent(crude, accuracy = 0.1)
  )

# 4) Crude mortality rate
crude_death_unique <- df_patient_flags %>%
  summarise(
    patients = n(),
    deaths   = sum(ever_death, na.rm = TRUE),
    crude    = deaths / patients,
    crude_pct = percent(crude, accuracy = 0.1)
  )

# 5) Distribution of mutually exclusive outcomes
tri_counts <- df_patient_flags %>%
  count(status_patient, name = "n") %>%
  mutate(pct = percent(n / sum(n), accuracy = 0.1))

print(crude_readmission_unique)
print(crude_death_unique)
print(tri_counts)


# =====================
# By infection types
# =====================
# 1) Patient-level collapse: keep infection type info
patient_type_flags <- df %>%
  mutate(recordid_base = stringr::str_remove(recordid, "_(VAP|BSI)$")) %>%
  group_by(recordid_base, infection_types) %>%
  summarise(
    ever_readm = as.integer(any(readm_event == 1, na.rm = TRUE)),
    .groups = "drop"
  )

# 2) Crude readmission by infection type
by_type_readm <- patient_type_flags %>%
  group_by(infection_types) %>%
  summarise(
    patients   = n_distinct(recordid_base),
    readmitted = sum(ever_readm, na.rm = TRUE),
    crude      = readmitted / patients,
    crude_pct  = percent(crude, accuracy = 0.1),
    .groups = "drop"
  ) %>%
  arrange(infection_types)

print(by_type_readm)


