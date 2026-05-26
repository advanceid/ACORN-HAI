# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr,
    dplyr,
    tidyr,
    AMR,
    tibble
  )
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
ast <- readRDS("data/clean_data/ast_all_index.RData")

# Add unique row id to preserve 1 row = 1 organism record
ast <- ast %>%
  mutate(.rowid = row_number())

# Build antibiotic-level AST dataset
ast_ab <- ast %>%
  select(-starts_with("ris_")) %>%
  rename(organism = org_names_all)

# Define metadata columns
meta_cols <- c(
  ".rowid",
  "recordid",
  "organism",
  "org",
  "org_oth",
  "f07m_bsic",
  "f07m_vapc",
  "ast_date",
  "spec_date",
  "f07m_hpd_adm_date_1",
  "f07m_hpd_adm_date_2",
  "redcap_repeat_instance",
  "org_names",
  "Class",
  "Methicillin",
  "org_new",
  "pathogen_group",
  "infection_types",
  "inf_onset",
  "age_new",
  "org_combined",
  "pathogen_group_combined",
  "pathogen_combined_types_ori",
  "pathogen_combined_types"
)

# Antibiotic columns
ab_cols <- setdiff(names(ast_ab), meta_cols)

# Convert AST coding:
# 1 = R, 2 = I, 3 = S, 5 = missing
ast_ab <- ast_ab %>%
  mutate(
    across(
      all_of(ab_cols),
      ~ case_when(
        . == 1 ~ "R",
        . == 2 ~ "I",
        . == 3 ~ "S",
        . == 5 ~ NA_character_,
        TRUE   ~ NA_character_
      )
    )
  ) %>%
  mutate(
    across(all_of(ab_cols), AMR::as.sir)
  )

# Standardise organism names
ast_ab <- ast_ab %>%
  mutate(
    mo = AMR::as.mo(organism)
  )

# Check unmatched organism names
unmatched_mo <- ast_ab %>%
  filter(is.na(mo) & !is.na(organism)) %>%
  distinct(organism)

# Keep rows with valid organism for downstream analyses
ast_ab_mdro <- ast_ab %>%
  filter(!is.na(organism), !is.na(mo))


# MDRO classification
mdro_result <- AMR::mdro(
  x = ast_ab_mdro,
  col_mo = "mo",
  guideline = "CMI 2012",
  pct_required_classes = 0.1,
  combine_SI = FALSE,
  verbose = FALSE
)

mdr_df <- ast_ab_mdro %>%
  transmute(
    .rowid,
    mdro = mdro_result,
    mdr = case_when(
      as.character(mdro_result) %in% c("MDR", "XDR", "PDR") ~ 1,
      as.character(mdro_result) == "Negative"               ~ 0,
      TRUE                                                  ~ NA_real_
    )
  )


# Drug-level AMR excluding intrinsic resistance
amr_long <- ast_ab_mdro %>%
  select(.rowid, mo, all_of(ab_cols)) %>%
  pivot_longer(
    cols = all_of(ab_cols),
    names_to = "ab",
    values_to = "sir"
  ) %>%
  mutate(
    intrinsic_r = AMR::mo_is_intrinsic_resistant(mo, ab),
    tested      = !is.na(sir),
    acquired_r  = sir == "R" & !intrinsic_r
  )

amr_ab_df <- amr_long %>%
  group_by(.rowid) %>%
  summarise(
    n_tested_ab_all          = sum(tested),
    n_tested_ab_nonintrinsic = sum(tested & !intrinsic_r, na.rm = TRUE),
    n_r_ab_all               = sum(sir == "R", na.rm = TRUE),
    n_r_ab_nonintrinsic      = sum(acquired_r, na.rm = TRUE),
    amr_ab = case_when(
      n_tested_ab_nonintrinsic == 0 ~ NA_real_,
      n_r_ab_nonintrinsic >= 1      ~ 1,
      TRUE                          ~ 0
    ),
    .groups = "drop"
  )


# Antibiotic -> class mapping
# ab_class() is deprecated, and amr_class() is not for this use.
# So map one antibiotic at a time.
get_one_class <- function(x) {
  out <- tryCatch(
    {
      info <- AMR::ab_info(x)
      if ("group" %in% names(info) && !is.na(info$group[1])) {
        as.character(info$group[1])
      } else if ("atc_group1" %in% names(info) && !is.na(info$atc_group1[1])) {
        as.character(info$atc_group1[1])
      } else if ("atc_group2" %in% names(info) && !is.na(info$atc_group2[1])) {
        as.character(info$atc_group2[1])
      } else {
        NA_character_
      }
    },
    error = function(e) NA_character_
  )
  out
}

ab_class_map <- tibble(
  ab = ab_cols,
  ab_class = vapply(ab_cols, get_one_class, character(1))
)

# Class-level AMR excluding intrinsic resistance
amr_class_long <- amr_long %>%
  left_join(ab_class_map, by = "ab") %>%
  filter(!is.na(ab_class))

class_level_df <- amr_class_long %>%
  group_by(.rowid, ab_class) %>%
  summarise(
    tested_class = any(tested & !intrinsic_r, na.rm = TRUE),
    r_class      = any(acquired_r, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(.rowid) %>%
  summarise(
    n_tested_class_nonintrinsic = sum(tested_class, na.rm = TRUE),
    n_r_class_nonintrinsic      = sum(r_class, na.rm = TRUE),
    amr_class = case_when(
      n_tested_class_nonintrinsic == 0 ~ NA_real_,
      n_r_class_nonintrinsic >= 1      ~ 1,
      TRUE                             ~ 0
    ),
    .groups = "drop"
  )


# Combine row-level outputs
ast_result <- ast_ab_mdro %>%
  select(.rowid) %>%
  left_join(mdr_df, by = ".rowid") %>%
  left_join(amr_ab_df, by = ".rowid") %>%
  left_join(class_level_df, by = ".rowid")


# Join back to original AST dataset
ast <- ast %>%
  left_join(ast_result, by = ".rowid") %>%
  select(-.rowid)


# Save
saveRDS(ast, "data/clean_data/ast_all_index.RData")
###
