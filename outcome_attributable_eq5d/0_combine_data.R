# Clear 
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr, 
    dplyr, 
    labelled,
    openxlsx,
    rlang
  )
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")

# Check for NA values in critical columns
columns_to_check <- c(
  "eq_5d_3l_1", "eq_5d_3l_2", "eq_5d_3l_3",
  "eq_5d_3l_4", "eq_5d_3l_5"
)

df_clean <- df[complete.cases(df[, columns_to_check]), ]

# Function to convert exposure to binary indicator
make_treated_indicator <- function(x) {
  if (is.numeric(x) || is.integer(x) || is.logical(x)) {
    return(as.integer(x == 1))
  }
  
  xv <- as.character(x)
  is_susc <- grepl("susceptible", xv, ignore.case = TRUE)
  is_resi <- grepl("resistant|car|van|meth|positive|yes", xv, ignore.case = TRUE)
  
  out <- integer(length(xv))
  out[is_resi & !is_susc] <- 1L
  
  return(out)
}

# Universal function to filter, clean data, and calculate weights
process_data_subset <- function(df, filter_column,
                                ps_covariates = c(
                                  "age_group_new", "sex", "country_region", 
                                  "country_income", "hpd_admreason", 
                                  "severity_score_scale", "icu_hd_ap",
                                  "pathogen_combined_types", "infection_types"
                                ),
                                weight_cap_quantile = 0.95,
                                prob_eps = 1e-6,
                                min_n_fit = 10) {
  
  if (!filter_column %in% colnames(df)) return(NULL)
  
  df <- df %>% 
    filter(!is.na(.data[[filter_column]])) %>%
    filter(age_new >= 12)
  
  keep_vars <- c(filter_column, ps_covariates)
  if (!all(keep_vars %in% names(df))) return(NULL)
  
  # Create treated indicator
  df$treated <- make_treated_indicator(df[[filter_column]])
  
  # Complete cases for propensity score model
  complete_idx <- complete.cases(df[, c("treated", ps_covariates), drop = FALSE])
  df_ps <- df[complete_idx, , drop = FALSE]
  
  # Skip if sample too small or no treatment variation
  if (nrow(df_ps) < min_n_fit || length(unique(df_ps$treated)) < 2L) return(NULL)
  
  # Fit propensity score model
  ps_formula <- as.formula(
    paste0("treated ~ ", paste(ps_covariates, collapse = " + "))
  )
  
  ps_mod <- tryCatch(
    glm(ps_formula, data = df_ps, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(ps_mod)) return(NULL)
  
  # Predict propensity score
  df$ps <- NA_real_
  df$ps[complete_idx] <- predict(ps_mod, newdata = df_ps, type = "response")
  
  # Avoid extreme 0 or 1 probabilities
  df$ps <- pmin(pmax(df$ps, prob_eps), 1 - prob_eps)
  
  # Marginal treatment probability
  p_treated <- mean(df_ps$treated == 1)
  
  # Stabilized IPTW weights
  df$wt <- NA_real_
  
  idx_t <- which(complete_idx & df$treated == 1)
  idx_c <- which(complete_idx & df$treated == 0)
  
  if (length(idx_t)) {
    df$wt[idx_t] <- p_treated / df$ps[idx_t]
  }
  
  if (length(idx_c)) {
    df$wt[idx_c] <- (1 - p_treated) / (1 - df$ps[idx_c])
  }
  
  # Truncate weights at 95th percentile
  cap <- quantile(df$wt, weight_cap_quantile, na.rm = TRUE)
  df$wt <- ifelse(!is.na(df$wt), pmin(df$wt, cap), NA_real_)
  
  # Relevel factors
  df$pathogen_combined_types <- droplevels(df$pathogen_combined_types)
  df$icu_hd_ap <- relevel(df$icu_hd_ap, ref = "No")
  
  # Set variable labels
  labels <- list(
    sex = "Sex",
    age_new = "Age (years)",
    age_group_new = "Age group",
    country_region = "Region",
    country_income = "World Bank income status",
    hpd_admreason = "Primary admission reason",
    icu_hd_ap = "Admission to ICU/HD at enrollment",
    pathogen_combined_types = "Pathogen type",
    infection_types = "Infection syndrome",
    comorbidities_CCI = "Charlson comorbidity index",
    sofa_score = "SOFA score",
    severity_score_scale = "Severity score of disease",
    fbis_score = "FBIS score",
    pitt_score = "PITT score",
    qpitt_score = "qPITT score",
    eq_5d_3l = "EQ-5D-3L score"
  )
  
  for (col in names(labels)) {
    if (col %in% colnames(df)) {
      df <- labelled::set_variable_labels(df, !!sym(col) := labels[[col]])
    }
  }
  
  # Keep rows with valid weights
  df <- df %>% 
    filter(!is.na(wt) & !is.na(ps)) %>%
    mutate(weight = wt)
  
  return(df)
}

# Function to relevel factors ending with "-susceptible"
relevel_factors_by_suffix <- function(df) {
  factor_columns <- colnames(df)[sapply(df, is.factor)] 
  
  for (col in factor_columns) {
    levels_col <- levels(df[[col]])
    
    if (any(grepl("-susceptible$", levels_col))) { 
      ref_level <- grep("-susceptible$", levels_col, value = TRUE)
      
      if (length(ref_level) == 1) {
        df[[col]] <- relevel(df[[col]], ref = ref_level)
      }
    }
  }
  
  return(df)
}

# Filter and clean subsets for different pathogens
pathogen_columns <- c(
  "aci_car", "ent_thir", "ent_car",
  "pse_car", "entc_van", "sa_meth"
)

df_used <- lapply(pathogen_columns, function(col) {
  sub_df <- process_data_subset(df_clean, col)
  
  if (!is.null(sub_df)) {
    sub_df <- relevel_factors_by_suffix(sub_df)
  }
  
  return(sub_df)
})

names(df_used) <- pathogen_columns

# All pathogens, subject more than 12 years old
df_clean2 <- df_clean %>% 
  filter(age_new >= 12)

# Save data
saveRDS(df_clean2, "data/att_eq_all.RData")
saveRDS(df_used, "data/att_eq.RData")
###
