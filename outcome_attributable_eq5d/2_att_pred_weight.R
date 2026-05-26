# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    magrittr
  )
})

# Working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_eq.RData")

# Pathogen/resistance variables
pathogen_columns <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth")

# Pretty labels
var_pretty <- c(
  aci_car   = "Carbapenem",
  ent_thir  = "Third-generation cephalosporin",
  ent_car   = "Carbapenem",
  pse_car   = "Carbapenem",
  entc_van  = "Vancomycin",
  sa_meth   = "Methicillin"
)

make_group_labels <- function(var_name) {
  drug <- var_pretty[[var_name]]
  if (is.null(drug) || is.na(drug)) drug <- var_name
  
  c(
    susceptible = paste0(drug, "-susceptible"),
    resistant   = paste0(drug, "-resistant")
  )
}

# Bootstrap marginal mean per observation with IPTW weights
bootstrap_mean_per_obs <- function(df_input, formula, var_name, by_var = NULL,
                                   n_boot = 500, seed = 123) {
  set.seed(seed)
  
  stopifnot(var_name %in% names(df_input))
  stopifnot(all(c("recordid", "eq_5d_3l", "weight") %in% names(df_input)))
  
  df_input <- as.data.frame(df_input) %>%
    filter(!is.na(weight), weight > 0)
  
  n_obs <- nrow(df_input)
  
  if (n_obs == 0) {
    stop("No observations with valid positive weights.")
  }
  
  pred_vec <- matrix(NA_real_, nrow = n_obs, ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    
    boot_idx <- sample.int(
      n_obs,
      size = n_obs,
      replace = TRUE
    )
    
    boot_sample <- df_input[boot_idx, , drop = FALSE]
    
    fit_b <- try(
      lm(formula, data = boot_sample, weights = weight),
      silent = TRUE
    )
    
    if (inherits(fit_b, "try-error")) next
    
    pred_b <- try(
      predict(fit_b, newdata = df_input),
      silent = TRUE
    )
    
    if (inherits(pred_b, "try-error")) next
    if (length(pred_b) != n_obs) next
    
    pred_vec[, b] <- as.numeric(pred_b)
  }
  
  valid_boots <- which(!is.na(pred_vec[1, ]))
  
  message(
    "Successful bootstraps for ", var_name,
    if (!is.null(by_var)) paste0(" [", by_var, "]") else "",
    ": ", length(valid_boots), "/", n_boot
  )
  
  if (length(valid_boots) < 20) {
    warning("Too few successful bootstrap replicates. Check model stability.")
  }
  
  out <- data.frame(
    recordid  = df_input$recordid,
    obs_id    = seq_len(n_obs),
    weight    = df_input$weight,
    mean_pred = rowMeans(pred_vec[, valid_boots, drop = FALSE], na.rm = TRUE),
    lower_ci  = apply(
      pred_vec[, valid_boots, drop = FALSE],
      1,
      stats::quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper_ci  = apply(
      pred_vec[, valid_boots, drop = FALSE],
      1,
      stats::quantile,
      probs = 0.975,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  out[[var_name]] <- df_input[[var_name]]
  
  if (!is.null(by_var) && by_var %in% names(df_input)) {
    out[[by_var]] <- df_input[[by_var]]
  }
  
  out
}

# Weighted marginal mean summary
mean_summary <- function(df_subset, group_label, subgroup_name, level_name) {
  
  if (is.null(df_subset) || nrow(df_subset) == 0) {
    return(NULL)
  }
  
  df_subset %>%
    summarise(
      adjusted_mean = weighted.mean(mean_pred, weight, na.rm = TRUE),
      lower_ci      = weighted.mean(lower_ci, weight, na.rm = TRUE),
      upper_ci      = weighted.mean(upper_ci, weight, na.rm = TRUE)
    ) %>%
    mutate(
      group    = group_label,
      subgroup = subgroup_name,
      level    = level_name
    ) %>%
    dplyr::select(subgroup, level, group, adjusted_mean, lower_ci, upper_ci)
}

# Count summary
count_summary <- function(df_result, df_mean, var_name) {
  
  if (is.null(df_mean) || nrow(df_mean) == 0) {
    return(NULL)
  }
  
  N_group <- df_result %>%
    group_by(.data[[var_name]]) %>%
    summarise(N = n_distinct(recordid), .groups = "drop")
  
  glab <- make_group_labels(var_name)
  
  N_group <- N_group %>%
    mutate(
      group = ifelse(
        grepl("sus", .data[[var_name]], ignore.case = TRUE),
        glab["susceptible"],
        glab["resistant"]
      )
    ) %>%
    dplyr::select(group, N)
  
  df_mean %>%
    left_join(N_group, by = "group") %>%
    dplyr::select(subgroup, level, group, N, adjusted_mean, lower_ci, upper_ci)
}

# Difference summary: resistant - susceptible
diff_summary <- function(df_mean, subgroup_name, level_name, var_name) {
  
  if (is.null(df_mean) || nrow(df_mean) == 0) {
    return(NULL)
  }
  
  glab <- make_group_labels(var_name)
  
  sus_row <- df_mean %>% filter(group == glab["susceptible"])
  res_row <- df_mean %>% filter(group == glab["resistant"])
  
  if (nrow(sus_row) == 0 || nrow(res_row) == 0) {
    return(NULL)
  }
  
  tibble(
    subgroup = subgroup_name,
    level    = level_name,
    pathogen = var_name,
    diff     = res_row$adjusted_mean - sus_row$adjusted_mean,
    lower_ci = res_row$lower_ci - sus_row$upper_ci,
    upper_ci = res_row$upper_ci - sus_row$lower_ci
  )
}

# Run one block
run_eq_block <- function(df_used, var_name, subgroup_var = NULL,
                         n_boot = 500, seed = 123) {
  
  glab <- make_group_labels(var_name)
  
  # Overall marginal model
  if (is.null(subgroup_var)) {
    
    fml <- as.formula(paste(
      "eq_5d_3l ~",
      var_name
    ))
    
    df_result <- bootstrap_mean_per_obs(
      df_input = df_used,
      formula  = fml,
      var_name = var_name,
      by_var   = NULL,
      n_boot   = n_boot,
      seed     = seed
    )
    
    df_mean <- bind_rows(
      mean_summary(
        df_result[grepl("sus", df_result[[var_name]], ignore.case = TRUE), ],
        glab["susceptible"],
        "overall",
        "overall"
      ),
      mean_summary(
        df_result[grepl("res", df_result[[var_name]], ignore.case = TRUE), ],
        glab["resistant"],
        "overall",
        "overall"
      )
    ) %>%
      mutate(pathogen = var_name, .before = 1)
    
    df_count <- count_summary(df_result, df_mean, var_name) %>%
      mutate(pathogen = var_name, .before = 1)
    
    df_diff <- diff_summary(df_mean, "overall", "overall", var_name)
    
    return(list(mean = df_mean, count = df_count, diff = df_diff))
  }
  
  # Subgroup marginal model with interaction
  fml <- as.formula(paste(
    "eq_5d_3l ~",
    subgroup_var, "*", var_name
  ))
  
  df_result <- bootstrap_mean_per_obs(
    df_input = df_used,
    formula  = fml,
    var_name = var_name,
    by_var   = subgroup_var,
    n_boot   = n_boot,
    seed     = seed
  )
  
  subgroup_levels <- unique(as.character(df_result[[subgroup_var]]))
  subgroup_levels <- subgroup_levels[!is.na(subgroup_levels)]
  
  mean_list  <- list()
  count_list <- list()
  diff_list  <- list()
  
  kk <- 1
  
  for (lv in subgroup_levels) {
    
    tmp <- df_result[as.character(df_result[[subgroup_var]]) == lv, , drop = FALSE]
    
    if (nrow(tmp) == 0) next
    
    df_mean_lv <- bind_rows(
      mean_summary(
        tmp[grepl("sus", tmp[[var_name]], ignore.case = TRUE), ],
        glab["susceptible"],
        subgroup_var,
        lv
      ),
      mean_summary(
        tmp[grepl("res", tmp[[var_name]], ignore.case = TRUE), ],
        glab["resistant"],
        subgroup_var,
        lv
      )
    ) %>%
      mutate(pathogen = var_name, .before = 1)
    
    df_count_lv <- count_summary(tmp, df_mean_lv, var_name) %>%
      mutate(pathogen = var_name, .before = 1)
    
    df_diff_lv <- diff_summary(df_mean_lv, subgroup_var, lv, var_name)
    
    mean_list[[kk]]  <- df_mean_lv
    count_list[[kk]] <- df_count_lv
    diff_list[[kk]]  <- df_diff_lv
    
    kk <- kk + 1
  }
  
  list(
    mean  = bind_rows(mean_list),
    count = bind_rows(count_list),
    diff  = bind_rows(diff_list)
  )
}

# Main loop
set.seed(123)

mean_results  <- list()
count_results <- list()
diff_results  <- list()

kk <- 1

for (i in seq_along(pathogen_columns)) {
  
  var_name <- pathogen_columns[i]
  df_used <- df_all[[i]]
  
  if (is.null(df_used) || nrow(df_used) == 0) {
    message("Skipping ", var_name, ": empty dataset.")
    next
  }
  
  df_used <- df_used %>%
    dplyr::select(where(~ !is.list(.))) %>%
    as.data.frame()
  
  needed_vars <- c(
    "recordid", "eq_5d_3l", "weight",
    "infection_types", "age_group_new", "country_income",
    var_name
  )
  
  needed_vars <- needed_vars[needed_vars %in% names(df_used)]
  df_used <- df_used[, needed_vars, drop = FALSE]
  
  if (!"weight" %in% names(df_used)) {
    message("Skipping ", var_name, ": weight variable not found.")
    next
  }
  
  key_vars <- c("recordid", "eq_5d_3l", "weight", var_name)
  key_vars <- key_vars[key_vars %in% names(df_used)]
  
  df_used <- df_used[complete.cases(df_used[, key_vars]), , drop = FALSE]
  df_used <- df_used[df_used$weight > 0, , drop = FALSE]
  
  if (nrow(df_used) == 0) {
    message("Skipping ", var_name, ": no complete cases with valid weights.")
    next
  }
  
  # Factors
  df_used[[var_name]] <- as.factor(df_used[[var_name]])
  
  fac_vars <- c("infection_types", "age_group_new", "country_income")
  fac_vars <- fac_vars[fac_vars %in% names(df_used)]
  
  for (v in fac_vars) {
    df_used[[v]] <- as.factor(df_used[[v]])
    df_used[[v]] <- droplevels(df_used[[v]])
  }
  
  df_used[[var_name]] <- droplevels(df_used[[var_name]])
  
  # Susceptible as reference
  lv <- levels(df_used[[var_name]])
  sus_level <- lv[grepl("sus", lv, ignore.case = TRUE)]
  
  if (length(sus_level) >= 1) {
    df_used[[var_name]] <- relevel(df_used[[var_name]], ref = sus_level[1])
  }
  
  # Overall
  out_overall <- run_eq_block(
    df_used = df_used,
    var_name = var_name,
    subgroup_var = NULL,
    n_boot = 500,
    seed = 123
  )
  
  mean_results[[kk]]  <- out_overall$mean
  count_results[[kk]] <- out_overall$count
  diff_results[[kk]]  <- out_overall$diff
  kk <- kk + 1
  
  # By infection syndrome
  if ("infection_types" %in% names(df_used)) {
    out_inf <- run_eq_block(
      df_used = df_used,
      var_name = var_name,
      subgroup_var = "infection_types",
      n_boot = 500,
      seed = 123
    )
    
    mean_results[[kk]]  <- out_inf$mean
    count_results[[kk]] <- out_inf$count
    diff_results[[kk]]  <- out_inf$diff
    kk <- kk + 1
  }
  
  # By age group
  if ("age_group_new" %in% names(df_used)) {
    out_age <- run_eq_block(
      df_used = df_used,
      var_name = var_name,
      subgroup_var = "age_group_new",
      n_boot = 500,
      seed = 123
    )
    
    mean_results[[kk]]  <- out_age$mean
    count_results[[kk]] <- out_age$count
    diff_results[[kk]]  <- out_age$diff
    kk <- kk + 1
  }
  
  # By country income
  if ("country_income" %in% names(df_used)) {
    out_inc <- run_eq_block(
      df_used = df_used,
      var_name = var_name,
      subgroup_var = "country_income",
      n_boot = 500,
      seed = 123
    )
    
    mean_results[[kk]]  <- out_inc$mean
    count_results[[kk]] <- out_inc$count
    diff_results[[kk]]  <- out_inc$diff
    kk <- kk + 1
  }
  
  message("Completed analysis for ", var_name)
}

# Combine all results
df_mean_all  <- bind_rows(mean_results)
df_count_all <- bind_rows(count_results)
df_diff_all  <- bind_rows(diff_results)

print(df_mean_all)
print(df_count_all)
print(df_diff_all)

# Save outputs
saveRDS(df_mean_all,  "data/eq5d3l_mean_all.RData")
saveRDS(df_count_all, "data/eq5d3l_count_all.RData")
saveRDS(df_diff_all,  "data/eq5d3l_diff_all.RData")
###