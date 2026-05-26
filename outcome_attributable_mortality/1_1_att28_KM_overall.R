# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 tidyr,
                 magrittr,
                 survival,
                 survminer,
                 ggplot2,
                 gridExtra,
                 ggpubr,
                 extrafont,
                 patchwork,
                 Cairo,
                 stringr,
                 grid)
})

# Set working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/att_first28.RData")

# Load fonts
loadfonts()

# Variables
var <- c("aci_car", "ent_thir", "ent_car", 
         "pse_car", "entc_van", "sa_meth", 
         "mdr_gnb", "mdr", "amr")

# Raw labels for group assignment
raw_labels <- list(
  c("CSA", "CRA"),
  c("3GCSE", "3GCRE"),
  c("CSE", "CRE"),
  c("CSP", "CRP"),
  c("VSE", "VRE"),
  c("MSSA", "MRSA"),
  c("Non-MDR-GNB", "MDR-GNB"),
  c("Non-MDR", "MDR"),
  c("Non-AMR", "AMR")
)

# Create padded labels for visual alignment
max_len <- max(nchar(unlist(raw_labels)))
labels <- lapply(raw_labels, function(x) str_pad(x, max_len, side = "left"))

# Initialize list to store plots
combined_plot <- list()

# Loop over each dataset
for (i in 1:9) {
  df <- df_all[[i]]
  ris_name <- var[i]
  
  # Define groupings
  df$ris_group <- factor(
    ifelse(grepl("resistant$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][2],
           ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][1], NA)),
    levels = raw_labels[[i]]
  )
  
  df$ris_label <- factor(
    ifelse(grepl("resistant$", df[[ris_name]], ignore.case = TRUE), labels[[i]][2],
           ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), labels[[i]][1], NA)),
    levels = labels[[i]]
  )
  
  df$time[df$time == 28 & df$event == 1] <- 27.999
  
  # Fit KM model
  fit <- survfit(Surv(time, event) ~ ris_group, data = df)
  
  # p-value
  p_model <- survdiff(Surv(time, event) ~ ris_group, data = df)
  p_value <- 1 - pchisq(p_model$chisq, df = length(p_model$n) - 1)
  p_label <- ifelse(
    p_value < 0.0001,
    "< 0.0001",
    paste0("= ", sprintf("%.4f", p_value))
  )
  
  # per-plot margins for y text/title 
  text_r  <- c(13, 15, 10, 10, 10, 10) 
  title_r <- c(15, 20, 15, 15, 15, 15) 
  text_r_table <- c(14, 12, 12, 12, 12, 12) 
  
  # Plot survival
  plt_fct <- ggsurvplot(
    fit,
    data = df,
    pval = FALSE,
    palette = c("#EFC000FF", "#0073c2FF"),
    conf.int = TRUE,
    linetype = 1,
    ggtheme = theme_minimal(base_family = "Times New Roman") +
      theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
        axis.text.x = element_text(color = "black", size = 16, family = "Times New Roman", hjust = 1.2),
        axis.text.y = element_text(color = "black", size = 16, family = "Times New Roman", 
                                   margin = margin(r = text_r[i])),
        axis.title.y = element_text(color = "black", size = 16, family = "Times New Roman", 
                                    margin = margin(r = title_r[i])),
        axis.title.x = element_text(color = "black", size = 16, family = "Times New Roman", margin = margin(t = 10)),
        panel.spacing = unit(0.6, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        strip.text.x = element_text(size = 16, margin = margin(t = 6, b = 6)),
        strip.background.x = element_rect(color = "black", fill = "gray90", linewidth = NA_real_)
      ),
    title = paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1]),
    xlab = "Follow-up time since infection onset (days)",
    ylab = "Survival probability",
    font.tickslab = c(16, "plain"),
    font.x = c(16, "plain"),
    font.y = c(16, "plain"),
    break.time.by = 7,
    xlim = c(0, 28),
    legend = "none"
  )
  
  # Risk table
  ggsurv <- ggsurvplot(fit, df, risk.table = TRUE, break.time.by = 7, xlim = c(0, 28), tables.theme = theme_minimal())
  
  tbl_data <- ggsurv$table$data %>%
    mutate(
      ris_group = gsub("ris_group=", "", strata),
      ris_group = factor(ris_group, levels = raw_labels[[i]])
    ) %>%
    filter(!is.na(n.risk), time <= 28)
  
  tbl_fct <- ggplot(tbl_data, aes(time, ris_group, color = ris_group)) +
    geom_text(aes(label = n.risk), size = 6, family = "Times New Roman") +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
      strip.background = element_rect(fill = "white", linewidth = 1),
      axis.text.y = element_text(color = "black", size = 16, family = "Times New Roman", 
                                 margin = margin(r = text_r_table[i])),
      axis.text.x = element_blank(),
      strip.text = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.6, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm")
    ) +
    scale_y_discrete(labels = labels[[i]]) +
    scale_x_continuous(
      limits = c(0, 28),
      breaks = c(0, 7, 14, 21, 28)
    ) +
    scale_color_manual(
      values = setNames(c("#EFC000FF", "#0073c2FF"), raw_labels[[i]]),
      labels = labels[[i]]
    ) +
    ggtitle(" ") + xlab("") + ylab("") + guides(color = "none")
  
  # Add title and p-value
  plt_fct$plot <- plt_fct$plot + ggtitle(paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1])) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, family = "Times New Roman")) +
    annotate("text", x = 5, y = 0.15,
             label = paste0("p ", p_label),
             size = 5.2, hjust = 1,
             family = "Times New Roman", fontface = "plain")
  
  # Combine plot + risk table
  grob_combined <- arrangeGrob(plt_fct$plot, tbl_fct, heights = c(7.5, 2.5))
  
  grob_padded <- gTree(
    children = gList(grob_combined),
    vp = viewport(layout = grid.layout(1, 1),
                  width = unit(1, "npc") - unit(1, "cm"),
                  height = unit(1, "npc") - unit(0.5, "cm"))
  )
  
  # Store the final plot in the list
  combined_plot[[i]] <- grob_padded

}

# Save
for (i in 1:6) {
  CairoPDF(file = paste0("output/sub_figure/att28_death_KM_0", i, "_overall.pdf"), width = 7, height = 5)
  grid.draw(combined_plot[[i]])
  dev.off()
}

# =========
# Helper: compute crude mortality by group
# =========
compute_crude <- function(df_all, var, raw_labels) {
  stopifnot(length(df_all) == length(var), length(var) == length(raw_labels))
  
  out_list <- vector("list", length(var))
  
  for (i in seq_along(var)) {
    df <- df_all[[i]]
    ris_name <- var[i]
    
    # Grouping consistent with your KM plots
    df$ris_group <- ifelse(
      grepl("resistant$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][2],
      ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][1], NA)
    )
    df$ris_group <- factor(df$ris_group, levels = raw_labels[[i]])
    
    # Summarize by ris_group
    tmp <- df %>%
      dplyr::filter(!is.na(ris_group)) %>%
      dplyr::group_by(ris_group) %>%
      dplyr::summarise(
        total   = dplyr::n(),
        deaths  = sum(event == 1, na.rm = TRUE),   # 28-day deaths
        crude   = deaths / total,
        crude_pct = scales::percent(crude, accuracy = 0.1),
        .groups = "drop"
      ) %>%
      dplyr::mutate(pattern = paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1]))
    
    out_list[[i]] <- tmp
  }
  
  dplyr::bind_rows(out_list) %>%
    dplyr::select(pattern, ris_group, total, deaths, crude, crude_pct)
}

crude_mortality <- compute_crude(df_all, var, raw_labels)
crude_mortality

# =========
# Helper: compute crude 28-day mortality by infection_types
# =========
compute_crude_by_infection_types <- function(df_all, var, raw_labels, inf_col = "infection_types") {
  stopifnot(length(df_all) == length(var), length(var) == length(raw_labels))
  
  out_list <- vector("list", length(var))
  
  for (i in seq_along(var)) {
    df <- df_all[[i]]
    ris_name <- var[i]
    
    if (!inf_col %in% names(df)) {
      stop(paste0("Column '", inf_col, "' not found in df_all[[", i, "]]. ",
                  "Available: ", paste(names(df), collapse = ", ")))
    }
    
    # Grouping consistent with KM plots
    df$ris_group <- ifelse(
      grepl("resistant$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][2],
      ifelse(grepl("susceptible$", df[[ris_name]], ignore.case = TRUE), raw_labels[[i]][1], NA)
    )
    df$ris_group <- factor(df$ris_group, levels = raw_labels[[i]])
    
    tmp <- df %>%
      dplyr::mutate(infection_types = as.character(.data[[inf_col]])) %>%
      dplyr::filter(!is.na(ris_group), !is.na(infection_types), infection_types != "") %>%
      dplyr::group_by(infection_types, ris_group) %>%
      dplyr::summarise(
        total     = dplyr::n(),
        deaths    = sum(event == 1, na.rm = TRUE),
        crude     = deaths / total,
        crude_pct = scales::percent(crude, accuracy = 0.1),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        pattern = paste(raw_labels[[i]][2], "versus", raw_labels[[i]][1]),
        ris_var = ris_name
      ) %>%
      dplyr::select(ris_var, pattern, infection_types, ris_group, total, deaths, crude, crude_pct)
    
    out_list[[i]] <- tmp
  }
  
  dplyr::bind_rows(out_list) %>%
    dplyr::arrange(ris_var, infection_types, ris_group)
}

crude_mortality_by_infection <- compute_crude_by_infection_types(df_all, var, raw_labels)

crude_mortality_by_infection

# =========================================================
# Make summary table (Overall + by infection_types) for i=1:6
# =========================================================

# 1) helper: format "xx.x% (d/t)"
fmt_cell <- function(deaths, total) {
  ifelse(
    is.na(total) | total == 0,
    NA_character_,
    paste0(sprintf("%.1f", 100 * deaths / total), "% (", deaths, "/", total, ")")
  )
}

# 2) clean infection types to avoid duplicates caused by spaces/case
clean_infection_type <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

# 3) build long summary for one organism index i
summarise_one <- function(df, ris_name, labels_pair, inf_col = "infection_types") {
  stopifnot(inf_col %in% names(df))
  
  df <- df %>%
    mutate(
      infection_types = clean_infection_type(.data[[inf_col]]),
      ris_group = case_when(
        grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[2],
        grepl("susceptible$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[1],
        TRUE ~ NA_character_
      ),
      ris_group = factor(ris_group, levels = labels_pair)
    ) %>%
    filter(!is.na(ris_group))
  
  # by infection_types
  by_inf <- df %>%
    filter(!is.na(infection_types)) %>%
    group_by(infection_types, ris_group) %>%
    summarise(
      total  = n(),
      deaths = sum(event == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cell = fmt_cell(deaths, total))
  
  # overall
  overall <- df %>%
    group_by(ris_group) %>%
    summarise(
      total  = n(),
      deaths = sum(event == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(infection_types = "Overall", cell = fmt_cell(deaths, total))
  
  bind_rows(overall, by_inf)
}

# 4) loop over first 6 organisms and bind
out_all6_long <- purrr::map_dfr(1:6, function(i) {
  df <- df_all[[i]]
  ris_name <- var[i]
  labs_i <- raw_labels[[i]]   # c(susceptible_label, resistant_label)
  
  summarise_one(df, ris_name, labs_i, inf_col = "infection_types") %>%
    mutate(
      pattern = paste(labs_i[2], "versus", labs_i[1]),
      ris_var = ris_name
    )
})

# 5) keep only the 4 blocks you want (optional, but usually matches your layout)
keep_blocks <- c("Overall", "VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")
out_all6_long <- out_all6_long %>%
  mutate(infection_types = factor(infection_types, levels = keep_blocks)) %>%
  filter(infection_types %in% keep_blocks)

# 6) wide table: two columns per pattern (susceptible vs resistant)
out_all6_wide <- out_all6_long %>%
  select(ris_var, pattern, infection_types, ris_group, cell) %>%
  tidyr::pivot_wider(
    names_from  = ris_group,
    values_from = cell,
    values_fn   = dplyr::first
  ) %>%
  arrange(ris_var, infection_types)

# 7) collapse to ONE column like: "Resistant ... vs. Susceptible ..."
#    raw_labels[[i]] = c(susceptible, resistant)
resistant_label_map <- purrr::map_chr(raw_labels[1:6], ~ .x[2])
susceptible_label_map <- purrr::map_chr(raw_labels[1:6], ~ .x[1])
names(resistant_label_map)   <- var[1:6]
names(susceptible_label_map) <- var[1:6]

out_all6_table <- out_all6_wide %>%
  mutate(
    Pattern = gsub("\\bversus\\b", "vs.", pattern),
    
    crude_mortality = dplyr::case_when(
      grepl("^CRA\\b", pattern)   ~ paste0(CRA,   " vs. ", CSA),
      grepl("^3GCRE\\b", pattern) ~ paste0(`3GCRE`, " vs. ", `3GCSE`),
      grepl("^CRE\\b", pattern)   ~ paste0(CRE,   " vs. ", CSE),
      grepl("^CRP\\b", pattern)   ~ paste0(CRP,   " vs. ", CSP),
      grepl("^VRE\\b", pattern)   ~ paste0(VRE,   " vs. ", VSE),
      grepl("^MRSA\\b", pattern)  ~ paste0(MRSA,  " vs. ", MSSA),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcomes = infection_types,
    Pattern,
    Crude_mortality = crude_mortality
  )

out_all6_table


# =========================================================
# Calculate crude mortality difference
# =========================================================

out_all6_table_diff <- out_all6_table %>%
  mutate(
    # extract two percentages
    pct_resistant   = as.numeric(stringr::str_extract(Crude_mortality, "^[0-9.]+")),
    pct_susceptible = as.numeric(stringr::str_extract(Crude_mortality, "(?<=vs\\. )[0-9.]+")),
    
    # difference
    Difference = pct_resistant - pct_susceptible,
    
    # format
    Difference = sprintf("%.1f%%", Difference)
  ) %>%
  select(-pct_resistant, -pct_susceptible)

out_all6_table_diff



# =========================================================
# Median days to death among 28-day deaths
# Groups: Non-AMR, AMR, MDR
# =========================================================

# helper: format median (IQR)
fmt_median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  
  paste0(
    sprintf(paste0("%.", digits, "f"), q[2]),
    " (IQR ",
    sprintf(paste0("%.", digits, "f"), q[1]),
    " to ",
    sprintf(paste0("%.", digits, "f"), q[3]),
    ")"
  )
}

# =========================================================
# Non-AMR and AMR from amr dataset (i = 9)
# =========================================================
df_amr_time <- df_all[[9]] %>%
  mutate(
    amr_group = case_when(
      grepl("resistant$", amr, ignore.case = TRUE)   ~ "AMR",
      grepl("susceptible$", amr, ignore.case = TRUE) ~ "Non-AMR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(amr_group),
    event == 1,
    !is.na(time)
  )


table(df_amr_time$amr_group)


amr_all <- df_amr_time[which(df_amr_time$amr_group == "AMR"),]
amr_death <- amr_all[which(amr_all$first28_death == 1),]

nrow(amr_death)
nrow(amr_all)

nrow(amr_death)/nrow(amr_all)



time_nonamr_amr <- df_amr_time %>%
  group_by(amr_group) %>%
  summarise(
    deaths = n(),
    median_days = median(time, na.rm = TRUE),
    q1 = quantile(time, 0.25, na.rm = TRUE),
    q3 = quantile(time, 0.75, na.rm = TRUE),
    median_iqr = fmt_median_iqr(time, digits = 1),
    .groups = "drop"
  )

print(time_nonamr_amr)

# =========================================================
# MDR from mdr dataset (i = 8)
# =========================================================
df_mdr_time <- df_all[[8]] %>%
  mutate(
    mdr_group = case_when(
      grepl("resistant$", mdr, ignore.case = TRUE)   ~ "MDR",
      grepl("susceptible$", mdr, ignore.case = TRUE) ~ "Non-MDR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    mdr_group == "MDR",
    event == 1,
    !is.na(time)
  )

time_mdr <- df_mdr_time %>%
  group_by(mdr_group) %>%
  summarise(
    deaths = n(),
    median_days = median(time, na.rm = TRUE),
    q1 = quantile(time, 0.25, na.rm = TRUE),
    q3 = quantile(time, 0.75, na.rm = TRUE),
    median_iqr = fmt_median_iqr(time, digits = 1),
    .groups = "drop"
  )

print(time_mdr)

# =========================================================
# Combine final summary
# =========================================================
time_summary_final <- bind_rows(
  time_nonamr_amr %>%
    transmute(group = amr_group, deaths, median_days, q1, q3, median_iqr),
  time_mdr %>%
    transmute(group = mdr_group, deaths, median_days, q1, q3, median_iqr)
) %>%
  mutate(
    group = factor(group, levels = c("Non-AMR", "AMR", "MDR"))
  ) %>%
  arrange(group)

print(time_summary_final)

# =========================================================
# Extract manuscript-ready text
# =========================================================
nonamr_time <- time_summary_final %>%
  filter(group == "Non-AMR") %>%
  pull(median_iqr)

amr_time <- time_summary_final %>%
  filter(group == "AMR") %>%
  pull(median_iqr)

mdr_time <- time_summary_final %>%
  filter(group == "MDR") %>%
  pull(median_iqr)


# =========================================================
# Crude 28-day mortality among MDR infections
# denominator = MDR patients within each subgroup
# =========================================================

df_mdr <- df_all[[8]] %>%
  mutate(
    mdr_group = case_when(
      grepl("resistant$", mdr, ignore.case = TRUE)   ~ "MDR",
      grepl("susceptible$", mdr, ignore.case = TRUE) ~ "Non-MDR",
      TRUE ~ NA_character_
    ),
    age_group_new = factor(
      age_group_new,
      levels = c("<1 year", "1–4 years", "5–14 years",
                 "15–49 years", "50–69 years", "≥70 years")
    )
  ) %>%
  filter(mdr_group == "MDR")

# =========================================================
# MDR crude mortality by 6 age groups
# numerator   = deaths among MDR patients in each subgroup
# denominator = MDR patients in each subgroup
# =========================================================
mdr_age6 <- df_mdr %>%
  filter(!is.na(age_group_new)) %>%
  group_by(age_group_new) %>%
  summarise(
    mdr_patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / mdr_patients,
    stat = paste0(sprintf("%.1f", pct), "% (", deaths, "/", mdr_patients, ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group_new)


all_age6 <- df_all[[9]] %>%
  filter(!is.na(age_group_new) & (amr == "-resistant")) %>%
  group_by(age_group_new) %>%
  summarise(
    patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / patients,
    stat = paste0(sprintf("%.1f", pct), "% (", deaths, "/", patients, ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group_new)


all_age6_country_income <- df_all[[9]] %>%
  filter(!is.na(country_income) & (amr == "-resistant")) %>%
  group_by(country_income) %>%
  summarise(
    patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / patients,
    stat = paste0(sprintf("%.1f", pct), "% (", deaths, "/", patients, ")"),
    .groups = "drop"
  ) %>%
  arrange(country_income)


# =========================================================
# Crude 28-day mortality by infection type, age, income
# AMR, MDR, CRA, 3GCRE, CRE
# =========================================================

get_crude_by_subgroup <- function(df, ris_name, resistant_label, subgroup_var) {
  
  df %>%
    mutate(
      ris_group = case_when(
        grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ resistant_label,
        TRUE ~ NA_character_
      ),
      subgroup_level = as.character(.data[[subgroup_var]])
    ) %>%
    filter(
      !is.na(ris_group),
      !is.na(subgroup_level),
      subgroup_level != ""
    ) %>%
    group_by(subgroup_level) %>%
    summarise(
      pathogen_label = resistant_label,
      total = n(),
      deaths = sum(event == 1, na.rm = TRUE),
      mortality = 100 * deaths / total,
      .groups = "drop"
    ) %>%
    mutate(
      subgroup = subgroup_var
    )
}

get_crude_all_subgroups <- function(df, ris_name, resistant_label) {
  
  bind_rows(
    get_crude_by_subgroup(df, ris_name, resistant_label, "infection_types"),
    get_crude_by_subgroup(df, ris_name, resistant_label, "age_group_new"),
    get_crude_by_subgroup(df, ris_name, resistant_label, "country_income")
  )
}

mort_subgroup_df <- bind_rows(
  get_crude_all_subgroups(df_all[[9]], "amr",      "AMR"),
  get_crude_all_subgroups(df_all[[8]], "mdr",      "MDR"),
  get_crude_all_subgroups(df_all[[1]], "aci_car",  "CRA"),
  get_crude_all_subgroups(df_all[[2]], "ent_thir", "3GCRE"),
  get_crude_all_subgroups(df_all[[3]], "ent_car",  "CRE")
)

# Add binomial CI
mort_subgroup_df <- mort_subgroup_df %>%
  rowwise() %>%
  mutate(
    ci = list(prop.test(deaths, total)$conf.int * 100),
    lo = ci[[1]],
    hi = ci[[2]],
    label = paste0(deaths, "/", total)
  ) %>%
  ungroup() %>%
  select(-ci)

# Clean labels
mort_subgroup_df <- mort_subgroup_df %>%
  mutate(
    subgroup_label = case_when(
      subgroup == "infection_types" ~ "Infection syndromes",
      subgroup == "age_group_new"   ~ "Age group",
      subgroup == "country_income"  ~ "World Bank income status"
    ),
    subgroup_level = case_when(
      subgroup_level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
      subgroup_level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
      subgroup_level == "Upper middle income" ~ "Upper middle\nincome",
      subgroup_level == "Lower middle income" ~ "Lower middle\nincome",
      TRUE ~ subgroup_level
    )
  )

# Factor order
mort_subgroup_df$subgroup_level <- factor(
  mort_subgroup_df$subgroup_level,
  levels = c(
    "VAP",
    "Hospital-\nacquired BSI",
    "Healthcare-\nassociated BSI",
    "<1 year", "1–4 years", "5–14 years",
    "15–49 years", "50–69 years", "≥70 years",
    "High income", "Upper middle\nincome", "Lower middle\nincome"
  )
)

mort_subgroup_df$subgroup_label <- factor(
  mort_subgroup_df$subgroup_label,
  levels = c(
    "Infection syndromes",
    "Age group",
    "World Bank income status"
  )
)

mort_subgroup_df$pathogen_label <- factor(
  mort_subgroup_df$pathogen_label,
  levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")
)

# Check data
mort_subgroup_df

















all_age3 <- df_all[[9]] %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group) %>%
  summarise(
    patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / patients,
    stat = paste0(sprintf("%.1f", pct), "% (", deaths, "/", patients, ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group)


all_age3 <- df_all[[9]] %>%
  filter(!is.na(age_new)) %>%
  mutate(
    age_group_new_new = case_when(
      age_new < 1 ~ "<1",
      age_new >= 1 & age_new <= 4 ~ "1–4",
      age_new >= 5 & age_new <= 17 ~ "5–17",
      age_new >= 18 & age_new <= 49 ~ "18–49",
      age_new >= 50 & age_new <= 69 ~ "50–69",
      age_new >= 70 ~ "≥70",
      TRUE ~ NA_character_
    ),
    age_group_new_new = factor(
      age_group_new_new,
      levels = c("<1", "1–4", "5–17", "18–49", "50–69", "≥70")
    )
  ) %>%
  group_by(age_group_new_new) %>%
  summarise(
    patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / patients,
    stat = paste0(sprintf("%.1f", pct), "% (", deaths, "/", patients, ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group_new_new)







# =========================================================
# MDR crude mortality by 3 age groups
# denominator = MDR patients within each subgroup
# =========================================================

mdr_age3 <- df_mdr %>%
  mutate(
    age_group_3 = case_when(
      age_group_new == "<1 year" ~ "<1 year",
      age_group_new %in% c("1–4 years", "5–14 years") ~ "1–14 years",
      age_group_new %in% c("15–49 years", "50–69 years", "≥70 years") ~ "≥15 years",
      TRUE ~ NA_character_
    ),
    age_group_3 = factor(age_group_3, levels = c("<1 year", "1–14 years", "≥15 years"))
  ) %>%
  filter(!is.na(age_group_3)) %>%
  group_by(age_group_3) %>%
  summarise(
    mdr_patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / mdr_patients,
    stat = paste0(
      sprintf("%.1f", pct),
      "% (", deaths, "/", mdr_patients, ")"
    ),
    .groups = "drop"
  ) %>%
  arrange(age_group_3)


# =========================================================
# MDR crude mortality by income level
# numerator   = deaths among MDR patients in each subgroup
# denominator = MDR patients in each subgroup
# =========================================================

mdr_income <- df_mdr %>%
  filter(!is.na(country_income)) %>%
  group_by(country_income) %>%
  summarise(
    mdr_patients = n(),
    deaths = sum(event == 1, na.rm = TRUE),
    pct = 100 * deaths / mdr_patients,
    stat = paste0(sprintf("%.1f", pct), "% (", deaths, "/", mdr_patients, ")"),
    .groups = "drop"
  ) %>%
  arrange(desc(pct))


# =========================================================
# Distribution of MDR infections by age group
# denominator = total number of MDR infections
# subgroup percentages sum to 100%
# =========================================================

df_mdr <- df_all[[8]] %>%
  mutate(
    mdr_group = case_when(
      grepl("resistant$", mdr, ignore.case = TRUE) ~ "MDR",
      grepl("susceptible$", mdr, ignore.case = TRUE) ~ "Non-MDR",
      TRUE ~ NA_character_
    ),
    age_group_new = factor(
      age_group_new,
      levels = c("<1 year", "1–4 years", "5–14 years",
                 "15–49 years", "50–69 years", "≥70 years")
    )
  )

# total MDR infections
total_mdr <- df_mdr %>%
  filter(mdr_group == "MDR") %>%
  nrow()

# by 6 age groups
mdr_age <- df_mdr %>%
  filter(mdr_group == "MDR", !is.na(age_group_new)) %>%
  group_by(age_group_new) %>%
  summarise(
    n = n(),
    pct = 100 * n / total_mdr,
    stat = paste0(sprintf("%.1f", pct), "% (", n, "/", total_mdr, ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group_new)

prop.test(mdr_age6$deaths, mdr_age6$mdr_patients)


# =========================================================
# Distribution of MDR infections by 3 age groups
# denominator = total number of MDR infections
# subgroup percentages sum to 100%
# =========================================================

mdr_age3 <- df_mdr %>%
  mutate(
    age_group_3 = case_when(
      age_group_new == "<1 year" ~ "<1 year",
      age_group_new %in% c("1–4 years", "5–14 years") ~ "1–14 years",
      age_group_new %in% c("15–49 years", "50–69 years", "≥70 years") ~ "≥15 years",
      TRUE ~ NA_character_
    ),
    age_group_3 = factor(age_group_3, levels = c("<1 year", "1–14 years", "≥15 years"))
  ) %>%
  filter(mdr_group == "MDR", !is.na(age_group_3)) %>%
  group_by(age_group_3) %>%
  summarise(
    n = n(),
    pct = 100 * n / total_mdr,
    stat = paste0(sprintf("%.1f", pct), "% (", n, "/", total_mdr, ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group_3)

# =========================================================
# Distribution of MDR infections by income level
# denominator = total number of MDR infections
# subgroup percentages sum to 100%
# =========================================================
mdr_income <- df_mdr %>%
  filter(mdr_group == "MDR", !is.na(country_income)) %>%
  group_by(country_income) %>%
  summarise(
    n = n(),
    pct = 100 * n / total_mdr,
    stat = paste0(sprintf("%.1f", pct), "% (", n, "/", total_mdr, ")"),
    .groups = "drop"
  ) %>%
  arrange(desc(n))
