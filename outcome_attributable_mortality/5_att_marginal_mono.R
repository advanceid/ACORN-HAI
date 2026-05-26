# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    tidyr,
    marginaleffects
  )
})

# Load data
pop <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
df_list <- readRDS("data/att_first28.RData")

df_list <- lapply(df_list, function(x) {
  x %>%
    filter(
      !is.na(pathogen_combined_types),
      pathogen_combined_types != "Polymicrobial"
    )
})

# resistance variables
pathogens <- c(
  "aci_car",
  "ent_thir",
  "ent_car",
  "pse_car",
  "entc_van",
  "sa_meth",
  "mdr_gnb",
  "mdr",
  "amr"
)

# helper: get susceptible / resistant counts
get_rs_counts <- function(data, var_name, subgroup = NULL) {
  
  keep_vars <- c(var_name, "first28_death")
  if (!is.null(subgroup)) keep_vars <- c(keep_vars, subgroup)
  
  data2 <- data %>%
    dplyr::select(all_of(keep_vars)) %>%
    mutate(
      exposure = case_when(
        grepl("resistant$", as.character(.data[[var_name]]), ignore.case = TRUE) ~ "Resistant",
        grepl("sus", as.character(.data[[var_name]]), ignore.case = TRUE) ~ "Susceptible",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(exposure))
  
  if (is.null(subgroup)) {
    
    out <- data2 %>%
      group_by(exposure) %>%
      summarise(
        n = n(),
        deaths = sum(first28_death == 1, na.rm = TRUE),
        risk = mean(first28_death, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = exposure,
        values_from = c(n, deaths, risk),
        names_sep = "_"
      ) %>%
      mutate(
        pathogen = var_name,
        subgroup = "overall",
        level = "overall"
      )
    
  } else {
    
    out <- data2 %>%
      group_by(.data[[subgroup]], exposure) %>%
      summarise(
        n = n(),
        deaths = sum(first28_death == 1, na.rm = TRUE),
        risk = mean(first28_death, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = exposure,
        values_from = c(n, deaths, risk),
        names_sep = "_"
      ) %>%
      mutate(
        pathogen = var_name,
        subgroup = subgroup,
        level = as.character(.data[[subgroup]])
      ) %>%
      select(pathogen, subgroup, level, everything(), -all_of(subgroup))
  }
  
  out
}

# Calculate AR from logistic models using df_list
results_list <- list()
k <- 1

for (i in seq_along(pathogens)) {
  
  var_name <- pathogens[i]
  df_used  <- df_list[[i]] %>%
    dplyr::select(where(~ !is.list(.)))
  
  # factor
  df_used[[var_name]] <- as.factor(df_used[[var_name]])
  
  if ("infection_types" %in% names(df_used)) {
    df_used$infection_types <- as.factor(df_used$infection_types)
  }
  if ("age_group_new" %in% names(df_used)) {
    df_used$age_group_new <- as.factor(df_used$age_group_new)
  }
  if ("country_income" %in% names(df_used)) {
    df_used$country_income <- as.factor(df_used$country_income)
  }
  
  # set susceptible as reference
  lv <- levels(df_used[[var_name]])
  sus_level <- lv[grepl("sus", lv, ignore.case = TRUE)]
  if (length(sus_level) == 1) {
    df_used[[var_name]] <- relevel(df_used[[var_name]], ref = sus_level)
  }
  
  # overall 
  fml_overall <- as.formula(
    paste0("first28_death ~ ", var_name)
  )
  
  model_overall <- glm(
    fml_overall,
    data = df_used,
    weights = wt,
    family = quasibinomial()
  )
  
  res_overall <- avg_comparisons(
    model_overall,
    variables = var_name,
    comparison = "difference",
    type = "response"
  )
  
  if (nrow(res_overall) > 0) {
    res_overall$pathogen <- var_name
    res_overall$subgroup <- "overall"
    res_overall$level <- "overall"
    res_overall$counts <- nrow(model.frame(model_overall))
    results_list[[k]] <- res_overall
    k <- k + 1
  }
  
  # infection_types
  if ("infection_types" %in% names(df_used)) {
    
    fml_inf <- as.formula(
      paste0("first28_death ~ infection_types * ", var_name)
    )
    
    model_inf <- glm(
      fml_inf,
      data = df_used,
      weights = wt,
      family = quasibinomial()
    )
    
    res_inf <- avg_comparisons(
      model_inf,
      variables = var_name,
      by = "infection_types",
      comparison = "difference",
      type = "response"
    )
    
    if (nrow(res_inf) > 0) {
      mf_inf <- model.frame(model_inf)
      
      n_inf <- mf_inf %>%
        dplyr::count(infection_types, name = "counts") %>%
        mutate(level = as.character(infection_types)) %>%
        select(level, counts)
      
      res_inf$pathogen <- var_name
      res_inf$subgroup <- "infection_types"
      res_inf$level <- as.character(res_inf$infection_types)
      res_inf <- res_inf %>% left_join(n_inf, by = "level")
      
      results_list[[k]] <- res_inf
      k <- k + 1
    }
  }
  
  # age_group_new
  if ("age_group_new" %in% names(df_used)) {
    
    fml_age <- as.formula(
      paste0("first28_death ~ age_group_new * ", var_name)
    )
    
    model_age <- glm(
      fml_age,
      data = df_used,
      weights = wt,
      family = quasibinomial()
    )
    
    res_age <- avg_comparisons(
      model_age,
      variables = var_name,
      by = "age_group_new",
      comparison = "difference",
      type = "response"
    )
    
    if (nrow(res_age) > 0) {
      mf_age <- model.frame(model_age)
      
      n_age <- mf_age %>%
        dplyr::count(age_group_new, name = "counts") %>%
        mutate(level = as.character(age_group_new)) %>%
        select(level, counts)
      
      res_age$pathogen <- var_name
      res_age$subgroup <- "age_group_new"
      res_age$level <- as.character(res_age$age_group_new)
      res_age <- res_age %>% left_join(n_age, by = "level")
      
      results_list[[k]] <- res_age
      k <- k + 1
    }
  }
  
  # country_income
  if ("country_income" %in% names(df_used)) {
    
    fml_inc <- as.formula(
      paste0("first28_death ~ country_income * ", var_name)
    )
    
    model_inc <- glm(
      fml_inc,
      data = df_used,
      weights = wt,
      family = quasibinomial()
    )
    
    res_inc <- avg_comparisons(
      model_inc,
      variables = var_name,
      by = "country_income",
      comparison = "difference",
      type = "response"
    )
    
    if (nrow(res_inc) > 0) {
      mf_inc <- model.frame(model_inc)
      
      n_inc <- mf_inc %>%
        dplyr::count(country_income, name = "counts") %>%
        mutate(level = as.character(country_income)) %>%
        select(level, counts)
      
      res_inc$pathogen <- var_name
      res_inc$subgroup <- "country_income"
      res_inc$level <- as.character(res_inc$country_income)
      res_inc <- res_inc %>% left_join(n_inc, by = "level")
      
      results_list[[k]] <- res_inc
      k <- k + 1
    }
  }
}

results_df <- bind_rows(results_list) %>%
  as.data.frame() %>%
  select(
    pathogen,
    subgroup,
    level,
    counts,
    contrast,
    estimate,
    conf.low,
    conf.high,
    statistic,
    p.value
  ) %>%
  mutate(
    AR = estimate,
    AR_lo = conf.low,
    AR_hi = conf.high,
    AR_CI = sprintf(
      "%.2f%% [%.2f, %.2f]",
      AR * 100, AR_lo * 100, AR_hi * 100
    )
  )


# Calculate subgroup-specific prevalence from pop
pre_list <- list()
m <- 1

for (var_name in pathogens) {
  
  pop[[var_name]] <- as.factor(pop[[var_name]])
  
  pop2 <- pop %>%
    mutate(
      resistant = grepl("-resistant$", as.character(.data[[var_name]]))
    )
  
  # overall
  pre_list[[m]] <- data.frame(
    pathogen   = var_name,
    subgroup   = "overall",
    level      = "overall",
    prevalence = mean(pop2$resistant, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  m <- m + 1
  
  # infection_types
  if ("infection_types" %in% names(pop2)) {
    tmp_inf <- pop2 %>%
      group_by(infection_types) %>%
      summarise(
        prevalence = mean(resistant, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        pathogen = var_name,
        subgroup = "infection_types",
        level = as.character(infection_types)
      ) %>%
      select(pathogen, subgroup, level, prevalence)
    
    pre_list[[m]] <- tmp_inf
    m <- m + 1
  }
  
  # age_group_new
  if ("age_group_new" %in% names(pop2)) {
    tmp_age <- pop2 %>%
      group_by(age_group_new) %>%
      summarise(
        prevalence = mean(resistant, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        pathogen = var_name,
        subgroup = "age_group_new",
        level = as.character(age_group_new)
      ) %>%
      select(pathogen, subgroup, level, prevalence)
    
    pre_list[[m]] <- tmp_age
    m <- m + 1
  }
  
  # country_income
  if ("country_income" %in% names(pop2)) {
    tmp_inc <- pop2 %>%
      group_by(country_income) %>%
      summarise(
        prevalence = mean(resistant, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        pathogen = var_name,
        subgroup = "country_income",
        level = as.character(country_income)
      ) %>%
      select(pathogen, subgroup, level, prevalence)
    
    pre_list[[m]] <- tmp_inc
    m <- m + 1
  }
}

pre_df <- bind_rows(pre_list)


# Calculate subgroup-specific Ipop from pop (same within subgroup across pathogens)
ipop_list <- list()

# overall
ipop_list[[1]] <- data.frame(
  subgroup = "overall",
  level    = "overall",
  Ipop     = mean(pop$first28_death, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# infection_types
if ("infection_types" %in% names(pop)) {
  ipop_list[[length(ipop_list) + 1]] <- pop %>%
    group_by(infection_types) %>%
    summarise(
      Ipop = mean(first28_death, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      subgroup = "infection_types",
      level = as.character(infection_types)
    ) %>%
    select(subgroup, level, Ipop)
}

# age_group_new
if ("age_group_new" %in% names(pop)) {
  ipop_list[[length(ipop_list) + 1]] <- pop %>%
    group_by(age_group_new) %>%
    summarise(
      Ipop = mean(first28_death, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      subgroup = "age_group_new",
      level = as.character(age_group_new)
    ) %>%
    select(subgroup, level, Ipop)
}

# country_income
if ("country_income" %in% names(pop)) {
  ipop_list[[length(ipop_list) + 1]] <- pop %>%
    group_by(country_income) %>%
    summarise(
      Ipop = mean(first28_death, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      subgroup = "country_income",
      level = as.character(country_income)
    ) %>%
    select(subgroup, level, Ipop)
}

ipop_df <- bind_rows(ipop_list)


# Diagnostic table: susceptible / resistant counts
diag_list <- list()
kk <- 1

for (i in seq_along(pathogens)) {
  
  var_name <- pathogens[i]
  df_used  <- df_list[[i]] %>%
    dplyr::select(where(~ !is.list(.)))
  
  # overall
  diag_list[[kk]] <- get_rs_counts(df_used, var_name, subgroup = NULL)
  kk <- kk + 1
  
  # infection_types
  if ("infection_types" %in% names(df_used)) {
    diag_list[[kk]] <- get_rs_counts(df_used, var_name, subgroup = "infection_types")
    kk <- kk + 1
  }
  
  # age_group_new
  if ("age_group_new" %in% names(df_used)) {
    diag_list[[kk]] <- get_rs_counts(df_used, var_name, subgroup = "age_group_new")
    kk <- kk + 1
  }
  
  # country_income
  if ("country_income" %in% names(df_used)) {
    diag_list[[kk]] <- get_rs_counts(df_used, var_name, subgroup = "country_income")
    kk <- kk + 1
  }
}

diag_df <- bind_rows(diag_list) %>%
  mutate(
    n_Resistant = dplyr::coalesce(n_Resistant, 0L),
    n_Susceptible = dplyr::coalesce(n_Susceptible, 0L),
    deaths_Resistant = dplyr::coalesce(deaths_Resistant, 0),
    deaths_Susceptible = dplyr::coalesce(deaths_Susceptible, 0),
    risk_Resistant = dplyr::coalesce(risk_Resistant, NA_real_),
    risk_Susceptible = dplyr::coalesce(risk_Susceptible, NA_real_)
  ) %>%
  mutate(
    resistant_rate = sprintf("%.2f%%", risk_Resistant * 100),
    susceptible_rate = sprintf("%.2f%%", risk_Susceptible * 100)
  )

# Merge and calculate PAF
final_df <- results_df %>%
  left_join(pre_df,  by = c("pathogen", "subgroup", "level")) %>%
  left_join(ipop_df, by = c("subgroup", "level")) %>%
  left_join(diag_df, by = c("pathogen", "subgroup", "level")) %>%
  # remove rows with N < 30
  filter(counts >= 30) %>%
  mutate(
    PAF = prevalence * AR / Ipop,
    PAF_lo = prevalence * AR_lo / Ipop,
    PAF_hi = prevalence * AR_hi / Ipop,
    
    num = paste0("N = ", counts),
    Prevalence = sprintf("%.2f%%", prevalence * 100),
    PAF_CI = sprintf(
      "%.2f%% [%.2f, %.2f]",
      PAF * 100, PAF_lo * 100, PAF_hi * 100
    )
  ) %>%
  select(
    pathogen,
    subgroup,
    level,
    num,
    Prevalence,
    n_Resistant,
    n_Susceptible,
    deaths_Resistant,
    deaths_Susceptible,
    AR_CI,
    PAF_CI
  )

print(final_df)

# Save
saveRDS(final_df, "data/att28_paf_marginal_mono.RData")
###

