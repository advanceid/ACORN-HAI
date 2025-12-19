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

# -----------------
# Check data
df_check <- df[which(df$first28_patient_days < 0), ] %>%
  select(recordid, inf_onset, mortality_date, first28_patient_days)

write.xlsx(df_check, "mortality_issues.xlsx")
# -----------------
# Delete first28_patient_days < 0 
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
  
  infection_types = "Infection syndromes",
  
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
         icu_hd_ap,
         infection_types, first28_death,
         aci_car, ent_thir, ent_car, 
         pse_car, entc_van, sa_meth) %>%
  tbl_strata(
    strata = first28_death,
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
  modify_table_styling(
    columns = label,
    rows = case_when(
      label == "Carbapenem-resistant Acinetobacter spp." ~ TRUE,
      label == "Third-generation cephalosporin-resistant Enterobacterales" ~ TRUE,
      label == "Carbapenem-resistant Enterobacterales" ~ TRUE,
      label == "Carbapenem-resistant Pseudomonas spp." ~ TRUE,
      label == "Vancomycin-resistant Enterococcus spp." ~ TRUE,
      label == "Methicillin-resistant Staphylococcus aureus" ~ TRUE,
      TRUE ~ FALSE
    ),
    footnote = "Absent indicated patients with susceptible pathogens or no pathogens for the specific category."
  ) %>%
  as_gt() %>%
  tab_options(footnotes.marks = c("*", "†", "#", "‡")) %>%
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
    1:5 ~ px(200),
    6:7 ~ px(70),
    8 ~ px(110),
    9 ~ px(130),
    10:12 ~ px(70),
    13 ~ px(110),
    14 ~ px(130)
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
gtsave(gt_table, filename = "output/table/baseline_table.html")

# ----------
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
# ===========================
# Clean record IDs and compute crude mortality
# ===========================

library(stringr)
library(scales)

# 1) Remove suffix (_VAP / _BSI) from recordid to get base patient ID
df_clean_id <- df %>%
  mutate(recordid_base = str_remove(recordid, "_(VAP|BSI)$"))

# 2) Collapse to unique patient-level dataset
#    Rule: if any episode has event == 1 (death), mark patient as dead
df_unique_id <- df_clean_id %>%
  group_by(recordid_base) %>%
  summarise(
    death = as.integer(any(event == 1, na.rm = TRUE)),
    .groups = "drop"
  )

# 3) Crude mortality at patient level
crude_mortality_unique <- df_unique_id %>%
  summarise(
    patients = n(),
    deaths   = sum(death, na.rm = TRUE),
    crude    = deaths / patients,
    crude_pct = percent(crude, accuracy = 0.1)
  )

print(crude_mortality_unique) 


# ===========================
# Crude mortality by infection type (episode level)
# ===========================
crude_by_infection <- df %>%
  group_by(infection_types) %>%
  summarise(
    episodes = n(),
    deaths   = sum(event == 1, na.rm = TRUE),
    crude    = deaths / episodes,
    crude_pct = percent(crude, accuracy = 0.1),
    .groups = "drop"
  ) %>%
  arrange(desc(crude))

print(crude_by_infection)

# Overall median time to death (among those who died)
death_time_unique <- df %>%
  dplyr::mutate(recordid_base = stringr::str_remove(recordid, "_(VAP|BSI)$")) %>%
  dplyr::filter(event == 1) %>%
  dplyr::group_by(recordid_base) %>%
  dplyr::summarise(
    time_to_death = min(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::summarise(
    n_patients = dplyr::n(),
    median_time = median(time_to_death, na.rm = TRUE),
    p25 = quantile(time_to_death, 0.25, na.rm = TRUE),
    p75 = quantile(time_to_death, 0.75, na.rm = TRUE)
  ) %>%
  dplyr::mutate(`Median [IQR] (days)` = sprintf("%.0f [%.0f, %.0f]", median_time, p25, p75)) %>%
  dplyr::select(n_patients, `Median [IQR] (days)`)

print(death_time_unique)

