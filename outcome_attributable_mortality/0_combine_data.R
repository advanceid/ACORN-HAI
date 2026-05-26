# Clear environment
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr,
                 dplyr,
                 labelled,
                 openxlsx,
                 rlang)
})

# Working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
ast <- readRDS("data/clean_data_RData/ast_all_index.RData")

# Keep only complete outcome cases
columns_to_check <- c("first28_patient_days", "first28_death")
df_clean <- df[complete.cases(df[, columns_to_check]), , drop = FALSE]

# Helper: relevel any factor ending with "susceptible"
relevel_factors_by_suffix <- function(dat) {
  if (is.null(dat) || !is.data.frame(dat)) return(dat)
  factor_cols <- names(dat)[vapply(dat, is.factor, logical(1))]
  for (col in factor_cols) {
    lvls <- levels(dat[[col]])
    if (!is.null(lvls) && any(grepl("susceptible$", lvls))) {
      ref <- grep("susceptible$", lvls, value = TRUE)
      if (length(ref) == 1) dat[[col]] <- stats::relevel(dat[[col]], ref = ref)
    }
  }
  dat
}

# =========================
# Helper: convert exposure to binary indicator
# Rule:
#   numeric/logical: 1 = treated
#   character/factor: "susceptible" = 0; resistant/car/van/meth/positive/yes = 1
# =========================
make_treated_indicator <- function(x) {
  if (is.numeric(x) || is.integer(x) || is.logical(x)) {
    return(as.integer(x == 1))
  }
  xv <- as.character(x)
  is_susc <- grepl("susceptible", xv, ignore.case = TRUE)
  is_resi <- grepl("resistant|car|van|meth|positive|yes", xv, ignore.case = TRUE)
  out <- integer(length(xv))
  out[is_resi & !is_susc] <- 1L
  out
}

# Core: preprocess each pathogen subset and compute IPTW weights
process_data_subset <- function(dat, filter_column,
                                ps_covariates = c("age_group_new","sex","country_region","country_income",
                                                  "hpd_admreason","severity_score_scale",
                                                  "icu_hd_ap","pathogen_combined_types","infection_types"),
                                weight_cap_quantile = 0.95,
                                prob_eps = 1e-6,
                                min_n_fit = 10) {
  if (is.null(dat) || !is.data.frame(dat)) return(NULL)
  if (!filter_column %in% names(dat))     return(NULL)
  
  # keep non-missing rows for this pathogen; clean follow-up time
  dat <- dat %>%
    dplyr::filter(!is.na(.data[[filter_column]])) %>%
    dplyr::mutate(first28_patient_days = ifelse(first28_patient_days <= 0, 1, first28_patient_days),
                  time  = first28_patient_days,
                  event = first28_death)
  
  keep_vars <- c(filter_column, ps_covariates)
  if (!all(keep_vars %in% names(dat))) return(NULL)
  
  # create treated indicator
  dat$treated <- make_treated_indicator(dat[[filter_column]])
  
  # only complete rows used for PS estimation
  complete_idx <- stats::complete.cases(dat[, c("treated", ps_covariates), drop = FALSE])
  dat_ps <- dat[complete_idx, , drop = FALSE]
  
  # skip if sample too small or no variation
  if (nrow(dat_ps) < min_n_fit || length(unique(dat_ps$treated)) < 2L) return(NULL)
  
  # fit PS model
  ps_formula <- stats::as.formula(
    paste0("treated ~ ", paste(ps_covariates, collapse = " + "))
  )
  ps_mod <- tryCatch(
    stats::glm(ps_formula, data = dat_ps, family = stats::binomial()),
    warning = function(w) w,
    error   = function(e) NULL
  )
  if (is.null(ps_mod) || inherits(ps_mod, "warning")) {
    return(NULL)
  }
  
  # predict and safely write back
  ps_vec <- stats::predict(ps_mod, newdata = dat_ps, type = "response")
  dat$ps <- NA_real_
  dat$ps[complete_idx] <- ps_vec
  
  # clip extreme probabilities
  dat$ps <- pmin(pmax(dat$ps, prob_eps), 1 - prob_eps)
  
  # marginal treatment probability
  p_treated <- mean(dat_ps$treated == 1)
  
  # stabilized IPTW weights
  dat$wt <- NA_real_
  idx_t <- which(complete_idx & dat$treated == 1)
  idx_c <- which(complete_idx & dat$treated == 0)
  if (length(idx_t)) dat$wt[idx_t] <- p_treated / dat$ps[idx_t]
  if (length(idx_c)) dat$wt[idx_c] <- (1 - p_treated) / (1 - dat$ps[idx_c])
  
  # truncate at 95th percentile
  cap <- stats::quantile(dat$wt, weight_cap_quantile, na.rm = TRUE)
  dat$wt <- ifelse(!is.na(dat$wt), pmin(dat$wt, cap), NA_real_)
  
  # factor releveling
  if ("icu_hd_ap" %in% names(dat) && is.factor(dat$icu_hd_ap)) {
    dat$icu_hd_ap <- stats::relevel(dat$icu_hd_ap, ref = "No")
  }
  if ("pathogen_combined_types" %in% names(dat) && is.factor(dat$pathogen_combined_types)) {
    dat$pathogen_combined_types <- droplevels(dat$pathogen_combined_types)
  }
  
  # labels (optional)
  labels <- list(
    sex = "Sex",
    age_new = "Age (years)",
    age_group_new = "Age group",
    country_region = "Region",
    country_income = "World Bank income status",
    hpd_admreason = "Primary admission reason",
    length_before_onset = "Days from admission to infection onset",
    icu_hd_ap = "ICU/HD admission at enrollment",
    pathogen_combined_types = "Pathogen type",
    comorbidities_CCI = "Charlson comorbidity index",
    sofa_score = "SOFA score",
    severity_score_scale = "Severity score",
    fbis_score = "FBIS score",
    pitt_score = "PITT score",
    qpitt_score = "qPITT score",
    eq_5d_3l = "EQ-5D-3L score",
    aci_car = "Acinetobacter spp.",
    ent_car = "Enterobacterales (carbapenem-resistant)",
    ent_thir = "Enterobacterales (3rd gen cephalosporin-resistant)",
    pse_car = "Pseudomonas spp.",
    entc_van = "Enterococcus spp.",
    sa_meth = "S. aureus",
    mdr_gnb = "MDR Gram-negative bacteria (≥3 classes)"
  )
  for (col in names(labels)) {
    if (col %in% names(dat)) dat <- labelled::set_variable_labels(dat, !!rlang::sym(col) := labels[[col]])
  }
  
  # return rows with valid ps and weights
  dat %>% dplyr::filter(!is.na(wt) & !is.na(ps))
}

# Apply to each pathogen column
pathogen_columns <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth", "mdr_gnb", "mdr", "amr")

df_used <- lapply(pathogen_columns, function(col) {
  sub_df <- process_data_subset(df_clean, col)
  if (!is.null(sub_df)) sub_df <- relevel_factors_by_suffix(sub_df)
  sub_df
})
names(df_used) <- pathogen_columns

# Save
saveRDS(df_used, "data/att_first28.RData")
###

