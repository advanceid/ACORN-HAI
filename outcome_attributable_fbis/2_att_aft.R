# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    MASS
  )
})

wd <- "./"
setwd(wd)

df_all <- readRDS("data/att_fbis.RData")

pathogen_columns <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth")

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

make_ordered_fbis <- function(x) {
  if (is.ordered(x)) return(x)
  x_chr <- as.character(x)
  lvls <- sort(unique(x_chr), na.last = TRUE)
  factor(x_chr, levels = lvls, ordered = TRUE)
}

# Bootstrap predicted probabilities per observation with IPTW weights
bootstrap_prob_per_obs <- function(df_input, formula, var_name, by_var = NULL,
                                   n_boot = 500, seed = 123) {
  set.seed(seed)
  
  stopifnot(var_name %in% names(df_input))
  stopifnot(all(c("recordid", "fbis_score", "weight") %in% names(df_input)))
  
  df_input <- as.data.frame(df_input) %>%
    filter(!is.na(weight), weight > 0)
  
  df_input$fbis_score <- make_ordered_fbis(df_input$fbis_score)
  
  n_obs <- nrow(df_input)
  score_levels <- levels(df_input$fbis_score)
  K <- length(score_levels)
  
  prob_array <- array(
    NA_real_,
    dim = c(n_obs, K, n_boot),
    dimnames = list(NULL, score_levels, NULL)
  )
  
  for (b in seq_len(n_boot)) {
    
    boot_idx <- sample.int(
      n_obs,
      size = n_obs,
      replace = TRUE
    )
    
    boot_sample <- df_input[boot_idx, , drop = FALSE]
    
    if (length(unique(boot_sample$fbis_score)) < K) next
    
    fit_b <- try(
      MASS::polr(
        formula,
        data = boot_sample,
        weights = weight,
        method = "logistic",
        Hess = TRUE
      ),
      silent = TRUE
    )
    
    if (inherits(fit_b, "try-error")) next
    
    pred_b <- try(
      predict(fit_b, newdata = df_input, type = "probs"),
      silent = TRUE
    )
    
    if (inherits(pred_b, "try-error")) next
    if (!all(score_levels %in% colnames(pred_b))) next
    
    pred_b <- pred_b[, score_levels, drop = FALSE]
    prob_array[, , b] <- as.matrix(pred_b)
  }
  
  valid_boots <- which(!is.na(prob_array[1, 1, ]))
  
  message(
    "Successful bootstraps for ", var_name,
    if (!is.null(by_var)) paste0(" [", by_var, "]") else "",
    ": ", length(valid_boots), "/", n_boot
  )
  
  if (length(valid_boots) < 20) {
    warning("Too few successful bootstrap replicates. Check model stability.")
  }
  
  out_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    probs_k <- prob_array[, k, valid_boots, drop = FALSE]
    
    df_k <- data.frame(
      recordid    = df_input$recordid,
      obs_id      = seq_len(n_obs),
      weight      = df_input$weight,
      score_level = score_levels[k],
      mean_prob   = rowMeans(probs_k, na.rm = TRUE),
      lower_ci    = apply(probs_k, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
      upper_ci    = apply(probs_k, 1, stats::quantile, probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    df_k[[var_name]] <- df_input[[var_name]]
    
    if (!is.null(by_var) && by_var %in% names(df_input)) {
      df_k[[by_var]] <- df_input[[by_var]]
    }
    
    out_list[[k]] <- df_k
  }
  
  bind_rows(out_list)
}

# Weighted marginal probability summary
prob_summary <- function(df_subset, group_label, subgroup_name, level_name) {
  df_subset %>%
    group_by(score_level) %>%
    summarise(
      mean_prob = weighted.mean(mean_prob, weight, na.rm = TRUE),
      lower_ci  = weighted.mean(lower_ci, weight, na.rm = TRUE),
      upper_ci  = weighted.mean(upper_ci, weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      fbis_score = score_level,
      group      = group_label,
      subgroup   = subgroup_name,
      level      = level_name
    ) %>%
    dplyr::select(subgroup, level, fbis_score, mean_prob, lower_ci, upper_ci, group)
}

count_summary <- function(df_result, df_plot, var_name) {
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
  
  df_plot %>%
    left_join(N_group, by = "group") %>%
    mutate(expected_count = round(mean_prob * N, 0)) %>%
    dplyr::select(subgroup, level, fbis_score, group, N, expected_count)
}

run_fbis_block <- function(df_used, var_name, subgroup_var = NULL,
                           n_boot = 500, seed = 123) {
  
  glab <- make_group_labels(var_name)
  
  # Overall marginal model
  if (is.null(subgroup_var)) {
    fml <- as.formula(paste(
      "fbis_score ~",
      var_name
    ))
    
    df_result <- bootstrap_prob_per_obs(
      df_input = df_used,
      formula  = fml,
      var_name = var_name,
      by_var   = NULL,
      n_boot   = n_boot,
      seed     = seed
    )
    
    df_plot <- bind_rows(
      prob_summary(
        df_result[grepl("sus", df_result[[var_name]], ignore.case = TRUE), ],
        glab["susceptible"],
        "overall",
        "overall"
      ),
      prob_summary(
        df_result[grepl("res", df_result[[var_name]], ignore.case = TRUE), ],
        glab["resistant"],
        "overall",
        "overall"
      )
    ) %>%
      mutate(pathogen = var_name)
    
    df_count <- count_summary(df_result, df_plot, var_name) %>%
      mutate(pathogen = var_name, .before = 1)
    
    return(list(plot = df_plot, count = df_count))
  }
  
  # Subgroup marginal model with interaction
  fml <- as.formula(paste(
    "fbis_score ~",
    subgroup_var, "*", var_name
  ))
  
  df_result <- bootstrap_prob_per_obs(
    df_input = df_used,
    formula  = fml,
    var_name = var_name,
    by_var   = subgroup_var,
    n_boot   = n_boot,
    seed     = seed
  )
  
  subgroup_levels <- unique(as.character(df_result[[subgroup_var]]))
  subgroup_levels <- subgroup_levels[!is.na(subgroup_levels)]
  
  df_plot_list  <- list()
  df_count_list <- list()
  kk <- 1
  
  for (lv in subgroup_levels) {
    tmp <- df_result[as.character(df_result[[subgroup_var]]) == lv, , drop = FALSE]
    if (nrow(tmp) == 0) next
    
    df_plot_lv <- bind_rows(
      prob_summary(
        tmp[grepl("sus", tmp[[var_name]], ignore.case = TRUE), ],
        glab["susceptible"],
        subgroup_var,
        lv
      ),
      prob_summary(
        tmp[grepl("res", tmp[[var_name]], ignore.case = TRUE), ],
        glab["resistant"],
        subgroup_var,
        lv
      )
    ) %>%
      mutate(pathogen = var_name)
    
    df_count_lv <- count_summary(tmp, df_plot_lv, var_name) %>%
      mutate(pathogen = var_name, .before = 1)
    
    df_plot_list[[kk]]  <- df_plot_lv
    df_count_list[[kk]] <- df_count_lv
    kk <- kk + 1
  }
  
  list(
    plot  = bind_rows(df_plot_list),
    count = bind_rows(df_count_list)
  )
}

# Main loop
set.seed(123)

plot_results  <- list()
count_results <- list()
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
    "recordid", "fbis_score", "weight",
    "infection_types", "age_group_new", "country_income",
    var_name
  )
  
  needed_vars <- needed_vars[needed_vars %in% names(df_used)]
  df_used <- df_used[, needed_vars, drop = FALSE]
  
  if (!"weight" %in% names(df_used)) {
    message("Skipping ", var_name, ": weight variable not found.")
    next
  }
  
  df_used <- df_used[
    complete.cases(df_used[, c("recordid", "fbis_score", "weight", var_name)]),
    ,
    drop = FALSE
  ]
  
  df_used <- df_used[df_used$weight > 0, , drop = FALSE]
  
  if (nrow(df_used) == 0) {
    message("Skipping ", var_name, ": no complete cases with valid weights.")
    next
  }
  
  df_used$fbis_score <- make_ordered_fbis(df_used$fbis_score)
  df_used[[var_name]] <- as.factor(df_used[[var_name]])
  
  fac_vars <- c("infection_types", "age_group_new", "country_income")
  fac_vars <- fac_vars[fac_vars %in% names(df_used)]
  
  for (v in fac_vars) {
    df_used[[v]] <- as.factor(df_used[[v]])
    df_used[[v]] <- droplevels(df_used[[v]])
  }
  
  df_used[[var_name]] <- droplevels(df_used[[var_name]])
  
  lv <- levels(df_used[[var_name]])
  sus_level <- lv[grepl("sus", lv, ignore.case = TRUE)]
  
  if (length(sus_level) >= 1) {
    df_used[[var_name]] <- relevel(df_used[[var_name]], ref = sus_level[1])
  }
  
  # Overall
  out_overall <- run_fbis_block(
    df_used = df_used,
    var_name = var_name,
    subgroup_var = NULL,
    n_boot = 500,
    seed = 123
  )
  
  plot_results[[kk]]  <- out_overall$plot
  count_results[[kk]] <- out_overall$count
  kk <- kk + 1
  
  # Infection syndrome
  if ("infection_types" %in% names(df_used)) {
    out_inf <- run_fbis_block(
      df_used = df_used,
      var_name = var_name,
      subgroup_var = "infection_types",
      n_boot = 500,
      seed = 123
    )
    
    plot_results[[kk]]  <- out_inf$plot
    count_results[[kk]] <- out_inf$count
    kk <- kk + 1
  }
  
  # Age group
  if ("age_group_new" %in% names(df_used)) {
    out_age <- run_fbis_block(
      df_used = df_used,
      var_name = var_name,
      subgroup_var = "age_group_new",
      n_boot = 500,
      seed = 123
    )
    
    plot_results[[kk]]  <- out_age$plot
    count_results[[kk]] <- out_age$count
    kk <- kk + 1
  }
  
  # Country income
  if ("country_income" %in% names(df_used)) {
    out_inc <- run_fbis_block(
      df_used = df_used,
      var_name = var_name,
      subgroup_var = "country_income",
      n_boot = 500,
      seed = 123
    )
    
    plot_results[[kk]]  <- out_inc$plot
    count_results[[kk]] <- out_inc$count
    kk <- kk + 1
  }
  
  message("Completed analysis for ", var_name)
}

df_plot_all  <- bind_rows(plot_results)
df_count_all <- bind_rows(count_results)

print(df_plot_all)
print(df_count_all)

saveRDS(df_plot_all,  "data/data_plot_all.RData")
saveRDS(df_count_all, "data/df_count_all.RData")


