# =========================================================
# FINAL COMBINED FIGURE
# (a) Crude 28-day mortality
# (b) Attributable mortality
# (c) Population attributable fraction
# =========================================================

rm(list = ls())

suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    tidyr,
    stringr,
    ggplot2,
    patchwork,
    Cairo,
    grid,
    survival
  )
})

# =========================================================
# Load data
# =========================================================

df_all <- readRDS("data/att_first28.RData")
final_df <- readRDS("data/att28_paf_marginal.RData")

# =========================================================
# COMMON THEME
# =========================================================

theme_my <- theme_minimal(
  base_size = 12,
  base_family = "Times"
) +
  theme(
    text = element_text(family = "Times", size = 12, color = "black"),
    plot.title = element_text(
      family = "Times",
      size = 12,
      face = "bold",
      color = "black",
      hjust = -0.07
    ),
    axis.title = element_text(family = "Times", size = 12, color = "black"),
    axis.text = element_text(family = "Times", size = 12, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length = unit(0.1, "cm"),
    strip.text = element_text(
      family = "Times",
      size = 12,
      face = "bold",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.25
    ),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.25),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

# =========================================================
# Colours
# =========================================================

my_cols <- c(
  "AMR"   = "#b85c5c",
  "MDR"   = "#d9a441",
  "CRA"   = "#6f8f6d",
  "3GCRE" = "#7d9a9b",
  "CRE"   = "#a28ab8"
)

my_fill <- c(
  "AMR"     = unname(my_cols["AMR"]),
  "Non-AMR" = unname(my_cols["AMR"]),
  "MDR"     = unname(my_cols["MDR"]),
  "Non-MDR" = unname(my_cols["MDR"]),
  "CRA"     = unname(my_cols["CRA"]),
  "CSA"     = unname(my_cols["CRA"]),
  "3GCRE"   = unname(my_cols["3GCRE"]),
  "3GCSE"   = unname(my_cols["3GCRE"]),
  "CRE"     = unname(my_cols["CRE"]),
  "CSE"     = unname(my_cols["CRE"])
)

# =========================================================
# PART A
# Crude mortality
# Resistant and susceptible circles
# =========================================================

get_pair_by_subgroup <- function(
    df,
    ris_name,
    labels_pair,
    pair_name,
    subgroup_var = NULL
) {
  
  if (is.null(subgroup_var)) {
    
    df2 <- df %>%
      mutate(
        group = case_when(
          grepl("susceptible$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[1],
          grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[2],
          TRUE ~ NA_character_
        ),
        level = "Overall",
        subgroup = "overall"
      ) %>%
      filter(!is.na(group))
    
  } else {
    
    df2 <- df %>%
      mutate(
        group = case_when(
          grepl("susceptible$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[1],
          grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[2],
          TRUE ~ NA_character_
        ),
        level = as.character(.data[[subgroup_var]]),
        subgroup = subgroup_var
      ) %>%
      filter(
        !is.na(group),
        !is.na(level),
        level != ""
      )
  }
  
  df2 %>%
    group_by(subgroup, level, group) %>%
    summarise(
      total = n(),
      deaths = sum(event == 1, na.rm = TRUE),
      mortality = 100 * deaths / total,
      .groups = "drop"
    ) %>%
    mutate(pair = pair_name)
}

get_pair_all_subgroups <- function(
    df,
    ris_name,
    labels_pair,
    pair_name
) {
  
  bind_rows(
    get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, NULL),
    get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "infection_types"),
    get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "age_group_new"),
    get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "country_income")
  )
}

mort_plot_df <- bind_rows(
  get_pair_all_subgroups(df_all[[9]], "amr", c("Non-AMR", "AMR"), "AMR"),
  get_pair_all_subgroups(df_all[[8]], "mdr", c("Non-MDR", "MDR"), "MDR"),
  get_pair_all_subgroups(df_all[[1]], "aci_car", c("CSA", "CRA"), "CRA"),
  get_pair_all_subgroups(df_all[[2]], "ent_thir", c("3GCSE", "3GCRE"), "3GCRE"),
  get_pair_all_subgroups(df_all[[3]], "ent_car", c("CSE", "CRE"), "CRE")
)

mort_plot_df <- mort_plot_df %>%
  mutate(
    subgroup_label = case_when(
      subgroup == "overall" ~ "Overall",
      subgroup == "infection_types" ~ "Infection syndromes",
      subgroup == "age_group_new" ~ "Age groups",
      subgroup == "country_income" ~ "World Bank income status"
    ),
    level = case_when(
      level == "overall" ~ "Overall",
      level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
      level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
      level == "Upper middle income" ~ "Upper middle\nincome",
      level == "Lower middle income" ~ "Lower middle\nincome",
      TRUE ~ level
    ),
    pair = factor(
      pair,
      levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")
    ),
    type = case_when(
      group %in% c(
        "Non-AMR",
        "Non-MDR",
        "CSA",
        "3GCSE",
        "CSE"
      ) ~ "Susceptible",
      TRUE ~ "Resistant"
    )
  )

mort_plot_df$level <- factor(
  mort_plot_df$level,
  levels = c(
    "Overall",
    "VAP",
    "Hospital-\nacquired BSI",
    "Healthcare-\nassociated BSI",
    "<1 year",
    "1–4 years",
    "5–14 years",
    "15–49 years",
    "50–69 years",
    "≥70 years",
    "High income",
    "Upper middle\nincome",
    "Lower middle\nincome"
  )
)

mort_plot_df$subgroup_label <- factor(
  mort_plot_df$subgroup_label,
  levels = c(
    "Overall",
    "Infection syndromes",
    "Age groups",
    "World Bank income status"
  )
)


p_mort <- ggplot(
  mort_plot_df,
  aes(
    x = level,
    y = mortality
  )
) +
  
  geom_line(
    aes(
      color = pair,
      group = interaction(pair, subgroup_label, level)
    ),
    position = position_dodge(width = 0.5),
    linewidth = 0.35,
    alpha = 0.7
  ) +
  
  geom_point(
    data = mort_plot_df %>% filter(type == "Susceptible"),
    aes(color = pair),
    position = position_dodge(width = 0.5),
    shape = 16,
    size = 1.8,
    alpha = 0.45
  ) +
  
  geom_point(
    data = mort_plot_df %>% filter(type == "Resistant"),
    aes(color = pair),
    position = position_dodge(width = 0.5),
    shape = 16,
    size = 1.8,
    alpha = 1
  ) +
  facet_grid(
    ~ subgroup_label,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_manual(
    values = my_cols,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Crude 28-day mortality\nacross s-r pathogen pairs (%)",
    color = NULL,
    alpha = NULL
  ) +
  theme_my +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, by = 20)
  )

# =========================================================
# PART B and C
# AR + PAF
# =========================================================

parse_ci <- function(x) {
  
  m <- stringr::str_match(
    x,
    "^\\s*([-]?[0-9.]+)%\\s*\\[\\s*([-]?[0-9.]+)\\s*,\\s*([-]?[0-9.]+)\\s*\\]\\s*$"
  )
  
  data.frame(
    est = as.numeric(m[, 2]),
    lo  = as.numeric(m[, 3]),
    hi  = as.numeric(m[, 4])
  )
}

if (!all(c("AR", "AR_lo", "AR_hi") %in% names(final_df))) {
  
  tmp <- parse_ci(final_df$AR_CI)
  
  final_df$AR    <- tmp$est
  final_df$AR_lo <- tmp$lo
  final_df$AR_hi <- tmp$hi
}

if (!all(c("PAF", "PAF_lo", "PAF_hi") %in% names(final_df))) {
  
  tmp <- parse_ci(final_df$PAF_CI)
  
  final_df$PAF    <- tmp$est
  final_df$PAF_lo <- tmp$lo
  final_df$PAF_hi <- tmp$hi
}

pathogen_map <- c(
  aci_car  = "CRA",
  ent_thir = "3GCRE",
  ent_car  = "CRE",
  pse_car  = "CRP",
  entc_van = "VRE",
  sa_meth  = "MRSA",
  mdr      = "MDR",
  amr      = "AMR"
)

final_df <- final_df %>%
  mutate(
    pathogen_label = recode(
      pathogen,
      !!!pathogen_map
    )
  )

plot_df <- final_df %>%
  filter(
    pathogen_label %in% c("AMR", "MDR", "CRA", "3GCRE", "CRE"),
    subgroup %in% c(
      "overall",
      "infection_types",
      "age_group_new",
      "country_income"
    ),
    !is.na(level)
  ) %>%
  mutate(
    subgroup_label = case_when(
      subgroup == "overall" ~ "Overall",
      subgroup == "infection_types" ~ "Infection syndromes",
      subgroup == "age_group_new" ~ "Age groups",
      subgroup == "country_income" ~ "World Bank income status"
    ),
    level = case_when(
      level == "overall" ~ "Overall",
      level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
      level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
      level == "Upper middle income" ~ "Upper middle\nincome",
      level == "Lower middle income" ~ "Lower middle\nincome",
      TRUE ~ level
    )
  )

plot_df$level <- factor(
  plot_df$level,
  levels = levels(mort_plot_df$level)
)

plot_df$subgroup_label <- factor(
  plot_df$subgroup_label,
  levels = levels(mort_plot_df$subgroup_label)
)

plot_df$pathogen_label <- factor(
  plot_df$pathogen_label,
  levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")
)

p_ar <- ggplot(
  plot_df,
  aes(
    x = level,
    y = AR,
    color = pathogen_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "grey80"
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = AR_lo,
      ymax = AR_hi
    ),
    position = position_dodge(width = 0.5),
    width = 0.25,
    linewidth = 0.3
  ) +
  facet_grid(
    ~ subgroup_label,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_manual(
    values = my_cols,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Attributable\nmortality (%)",
    color = NULL
  ) +
  theme_my +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  scale_y_continuous(
    breaks = seq(-60, 60, by = 20)
  )

p_paf <- ggplot(
  plot_df,
  aes(
    x = level,
    y = PAF,
    color = pathogen_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "grey80"
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = PAF_lo,
      ymax = PAF_hi
    ),
    position = position_dodge(width = 0.5),
    width = 0.25,
    linewidth = 0.3
  ) +
  facet_grid(
    ~ subgroup_label,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_manual(
    values = my_cols,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  ) +
  labs(
    x = NULL,
    y = "Population attributable\nfraction (%)",
    color = NULL
  ) +
  theme_my +
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Times",
      size = 10,
      color = "black",
      margin = margin(l = 4, r = 10)
    ),
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.box.margin = margin(t = -5),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(
      family = "Times",
      size = 12,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      color = "black",
      margin = margin(t = 3)
    )
  ) +
  scale_y_continuous(
    breaks = seq(-40, 100, by = 20)
  )


# =========================================================
# CUSTOM LEGEND
# =========================================================

legend_df <- data.frame(
  label = c(
    "AMR", "MDR", "CRA", "3GCRE", "CRE",
    "Non-AMR", "Non-MDR", "CSA", "3GCSE", "CSE"
  ),
  pair = c(
    "AMR", "MDR", "CRA", "3GCRE", "CRE",
    "AMR", "MDR", "CRA", "3GCRE", "CRE"
  ),
  type = c(
    rep("Resistant", 5),
    rep("Susceptible", 5)
  ),
  x = c(
    3, 3.3, 3.6, 3.8, 4.05,
    3, 3.3, 3.6, 3.8, 4.05
  ),
  y = c(rep(2, 5), rep(1, 5))
)

p_legend <- ggplot(
  legend_df,
  aes(x = x, y = y)
) +
  
  geom_segment(
    aes(
      x = x - 0.035,
      xend = x + 0.035,
      y = y,
      yend = y,
      color = pair,
      alpha = type
    ),
    linewidth = 0.3
  ) +
  
  geom_point(
    aes(
      color = pair,
      alpha = type
    ),
    shape = 16,
    size = 2.2
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    color = "black",
    nudge_x = 0.05,
    hjust = 0,
    family = "Times",
    size = 3.2
  ) +
  
  scale_color_manual(
    values = my_cols
  ) +
  
  scale_alpha_manual(
    values = c(
      "Resistant" = 1,
      "Susceptible" = 0.35
    )
  ) +
  
  coord_cartesian(
    xlim = c(2, 5),
    ylim = c(0.6, 2.4),
    clip = "off"
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "none",
    plot.margin = margin(
      t = -5,
      r = 10,
      b = 0,
      l = 10
    )
  )

# =========================================================
# COMBINE
# =========================================================

final_plot <- (p_mort / p_ar / p_paf / p_legend) +
  plot_layout(
    heights = c(1, 1, 1, 0.18)
  )


CairoPDF(
  file = "output/figure 3.pdf",
  width = 15,
  height = 8
)

print(final_plot)

dev.off()






# # =========================================================
# # FINAL COMBINED FIGURE
# # A: Crude mortality as soft circles
# # B: Attributable mortality
# # C: Population attributable fraction
# # =========================================================
# 
# rm(list = ls())
# 
# suppressPackageStartupMessages({
#   require(pacman)
#   pacman::p_load(
#     dplyr, tidyr, stringr, ggplot2, patchwork,
#     Cairo, grid, survival, scales
#   )
# })
# 
# df_all <- readRDS("data/att_first28.RData")
# final_df <- readRDS("data/att28_paf_marginal.RData")
# 
# # =========================================================
# # COLOURS
# # =========================================================
# 
# my_cols <- c(
#   "AMR"   = "#b85c5c",
#   "MDR"   = "#d9a441",
#   "CRA"   = "#6f8f6d",
#   "3GCRE" = "#7d9a9b",
#   "CRE"   = "#a28ab8"
# )
# 
# my_fill <- c(
#   "AMR"     = unname(my_cols["AMR"]),
#   "Non-AMR" = unname(my_cols["AMR"]),
#   "MDR"     = unname(my_cols["MDR"]),
#   "Non-MDR" = unname(my_cols["MDR"]),
#   "CRA"     = unname(my_cols["CRA"]),
#   "CSA"     = unname(my_cols["CRA"]),
#   "3GCRE"   = unname(my_cols["3GCRE"]),
#   "3GCSE"   = unname(my_cols["3GCRE"]),
#   "CRE"     = unname(my_cols["CRE"]),
#   "CSE"     = unname(my_cols["CRE"])
# )
# 
# # =========================================================
# # THEME
# # =========================================================
# 
# theme_my <- theme_minimal(base_size = 12, base_family = "Times") +
#   theme(
#     text = element_text(family = "Times", size = 12, color = "black"),
#     plot.title = element_text(family = "Times", size = 14, face = "bold", color = "black"),
#     axis.title = element_text(family = "Times", size = 12, color = "black"),
#     axis.text = element_text(family = "Times", size = 12, color = "black"),
#     axis.ticks = element_line(color = "black", linewidth = 0.3),
#     axis.ticks.length = unit(-0.1, "cm"),
#     strip.text = element_text(family = "Times", size = 12, color = "black"),
#     strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.35),
#     panel.grid.major = element_line(color = "grey94", linewidth = 0.22),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
#     legend.background = element_rect(fill = "white", color = NA),
#     legend.key = element_rect(fill = "white", color = NA)
#   )
# 
# # =========================================================
# # PART A DATA
# # =========================================================
# 
# get_pair_by_subgroup <- function(df, ris_name, labels_pair, pair_name, subgroup_var) {
#   df %>%
#     mutate(
#       group = case_when(
#         grepl("susceptible$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[1],
#         grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[2],
#         TRUE ~ NA_character_
#       ),
#       level = as.character(.data[[subgroup_var]]),
#       subgroup = subgroup_var
#     ) %>%
#     filter(!is.na(group), !is.na(level), level != "") %>%
#     group_by(subgroup, level, group) %>%
#     summarise(
#       total = n(),
#       deaths = sum(event == 1, na.rm = TRUE),
#       mortality = 100 * deaths / total,
#       .groups = "drop"
#     ) %>%
#     mutate(pair = pair_name)
# }
# 
# get_pair_all_subgroups_A <- function(df, ris_name, labels_pair, pair_name) {
#   bind_rows(
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "infection_types"),
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "country_income"),
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "age_group_new")
#   )
# }
# 
# mort_plot_df <- bind_rows(
#   get_pair_all_subgroups_A(df_all[[9]], "amr", c("Non-AMR", "AMR"), "AMR"),
#   get_pair_all_subgroups_A(df_all[[8]], "mdr", c("Non-MDR", "MDR"), "MDR"),
#   get_pair_all_subgroups_A(df_all[[1]], "aci_car", c("CSA", "CRA"), "CRA"),
#   get_pair_all_subgroups_A(df_all[[2]], "ent_thir", c("3GCSE", "3GCRE"), "3GCRE"),
#   get_pair_all_subgroups_A(df_all[[3]], "ent_car", c("CSE", "CRE"), "CRE")
# ) %>%
#   mutate(
#     subgroup_label = case_when(
#       subgroup == "infection_types" ~ "Infection syndromes",
#       subgroup == "country_income" ~ "World Bank income status",
#       subgroup == "age_group_new" ~ "Age groups"
#     ),
#     level = case_when(
#       level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
#       level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
#       level == "Upper middle income" ~ "Upper middle\nincome",
#       level == "Lower middle income" ~ "Lower middle\nincome",
#       TRUE ~ level
#     ),
#     pair = factor(pair, levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")),
#     group = factor(
#       group,
#       levels = c(
#         "AMR", "Non-AMR",
#         "MDR", "Non-MDR",
#         "CRA", "CSA",
#         "3GCRE", "3GCSE",
#         "CRE", "CSE"
#       )
#     ),
#     alpha_group = case_when(
#       group %in% c("Non-AMR", "Non-MDR", "CSA", "3GCSE", "CSE") ~ 0.42,
#       TRUE ~ 0.68
#     )
#   )
# 
# mort_plot_df$level <- factor(
#   mort_plot_df$level,
#   levels = c(
#     "VAP",
#     "Hospital-\nacquired BSI",
#     "Healthcare-\nassociated BSI",
#     "High income",
#     "Upper middle\nincome",
#     "Lower middle\nincome",
#     "<1 year",
#     "1–4 years",
#     "5–14 years",
#     "15–49 years",
#     "50–69 years",
#     "≥70 years"
#   )
# )
# 
# # =========================================================
# # PART A PLOT
# # =========================================================
# 
# make_block_title <- function(title) {
#   ggplot() +
#     annotate(
#       "rect",
#       xmin = -Inf, xmax = Inf,
#       ymin = -Inf, ymax = Inf,
#       fill = "grey90",
#       color = "black",
#       linewidth = 0.35
#     ) +
#     annotate(
#       "text",
#       x = 0.5, y = 0.5,
#       label = title,
#       family = "Times",
#       size = 4
#     ) +
#     xlim(0, 1) +
#     ylim(0, 1) +
#     theme_void() +
#     theme(plot.margin = margin(t = 0, r = 0, b = -5, l = 0))
# }
# 
# make_panel_title <- function(title) {
#   ggplot() +
#     annotate(
#       "text",
#       x = -0.12, y = 0.5,
#       label = title,
#       family = "Times",
#       fontface = "bold",
#       size = 4.2,
#       hjust = 0
#     ) +
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
#     theme_void() +
#     theme(plot.margin = margin(t = 6, r = 0, b = 1, l = 0))
# }
# 
# plot_A_block <- function(dat, show_y = FALSE) {
#   
#   p <- ggplot(
#     dat,
#     aes(
#       x = pair,
#       y = mortality,
#       fill = group,
#       alpha = alpha_group
#     )
#   ) +
#     geom_point(
#       position = position_dodge(width = 0.22),
#       shape = 21,
#       size = 2.6,
#       stroke = 0.28,
#       color = "grey55"
#     ) +
#     scale_alpha_identity() +
#     facet_grid(~ level, scales = "free_x", space = "free_x") +
#     scale_fill_manual(values = my_fill, drop = FALSE) +
#     scale_y_continuous(
#       limits = c(0, 90),
#       breaks = seq(0, 80, 20),
#       labels = function(x) paste0(x, "%"),
#       expand = expansion(mult = c(0, 0.04))
#     ) +
#     labs(x = NULL, fill = NULL) +
#     theme_my +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       strip.text.x = element_text(family = "Times", size = 12, color = "black"),
#       panel.spacing.x = unit(0.03, "cm"),
#       plot.margin = margin(t = 0, r = 2, b = 0, l = 2)
#     )
#   
#   if (show_y) {
#     p <- p +
#       labs(y = "Crude 28-day\nmortality (%)") +
#       theme(
#         axis.title.y = element_text(family = "Times", size = 12, color = "black"),
#         axis.text.y = element_text(family = "Times", size = 12, color = "black"),
#         axis.ticks.y = element_line(color = "black", linewidth = 0.25)
#       )
#   } else {
#     p <- p +
#       labs(y = NULL) +
#       theme(
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#       )
#   }
#   
#   p
# }
# 
# p_A_syndrome_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "Infection syndromes"),
#   show_y = TRUE
# )
# 
# p_A_income_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "World Bank income status"),
#   show_y = FALSE
# )
# 
# p_A_age_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "Age groups"),
#   show_y = TRUE
# )
# 
# p_A_syndrome <- patchwork::wrap_plots(
#   make_block_title("Infection syndromes"),
#   p_A_syndrome_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_A_income <- patchwork::wrap_plots(
#   make_block_title("World Bank income status"),
#   p_A_income_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_A_age <- patchwork::wrap_plots(
#   make_block_title("Age groups"),
#   p_A_age_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_mort_body <- (p_A_syndrome | p_A_income) /
#   patchwork::plot_spacer() /
#   p_A_age +
#   plot_layout(heights = c(1, 0.002, 1.05))
# 
# p_mort <- patchwork::wrap_plots(
#   make_panel_title("(a) Crude mortality by resistant–susceptible pairs"),
#   p_mort_body,
#   ncol = 1,
#   heights = c(0.08, 1)
# )
# 
# # =========================================================
# # PART B/C DATA
# # =========================================================
# 
# parse_ci <- function(x) {
#   m <- stringr::str_match(
#     x,
#     "^\\s*([-]?[0-9.]+)%\\s*\\[\\s*([-]?[0-9.]+)\\s*,\\s*([-]?[0-9.]+)\\s*\\]\\s*$"
#   )
#   data.frame(
#     est = as.numeric(m[, 2]),
#     lo  = as.numeric(m[, 3]),
#     hi  = as.numeric(m[, 4])
#   )
# }
# 
# if (!all(c("AR", "AR_lo", "AR_hi") %in% names(final_df))) {
#   tmp <- parse_ci(final_df$AR_CI)
#   final_df$AR <- tmp$est
#   final_df$AR_lo <- tmp$lo
#   final_df$AR_hi <- tmp$hi
# }
# 
# if (!all(c("PAF", "PAF_lo", "PAF_hi") %in% names(final_df))) {
#   tmp <- parse_ci(final_df$PAF_CI)
#   final_df$PAF <- tmp$est
#   final_df$PAF_lo <- tmp$lo
#   final_df$PAF_hi <- tmp$hi
# }
# 
# pathogen_map <- c(
#   aci_car = "CRA",
#   ent_thir = "3GCRE",
#   ent_car = "CRE",
#   pse_car = "CRP",
#   entc_van = "VRE",
#   sa_meth = "MRSA",
#   mdr = "MDR",
#   amr = "AMR"
# )
# 
# plot_df <- final_df %>%
#   mutate(pathogen_label = recode(pathogen, !!!pathogen_map)) %>%
#   filter(
#     pathogen_label %in% c("AMR", "MDR", "CRA", "3GCRE", "CRE"),
#     subgroup %in% c("infection_types", "age_group_new", "country_income"),
#     !is.na(level)
#   ) %>%
#   mutate(
#     subgroup_label = case_when(
#       subgroup == "infection_types" ~ "Infection syndromes",
#       subgroup == "age_group_new" ~ "Age groups",
#       subgroup == "country_income" ~ "World Bank income status"
#     ),
#     level = case_when(
#       level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
#       level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
#       level == "Upper middle income" ~ "Upper middle\nincome",
#       level == "Lower middle income" ~ "Lower middle\nincome",
#       TRUE ~ level
#     ),
#     level = factor(
#       level,
#       levels = c(
#         "VAP",
#         "Hospital-\nacquired BSI",
#         "Healthcare-\nassociated BSI",
#         "<1 year",
#         "1–4 years",
#         "5–14 years",
#         "15–49 years",
#         "50–69 years",
#         "≥70 years",
#         "High income",
#         "Upper middle\nincome",
#         "Lower middle\nincome"
#       )
#     ),
#     subgroup_label = factor(
#       subgroup_label,
#       levels = c(
#         "Infection syndromes",
#         "Age groups",
#         "World Bank income status"
#       )
#     ),
#     pathogen_label = factor(
#       pathogen_label,
#       levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")
#     )
#   )
# 
# # =========================================================
# # PART B
# # =========================================================
# 
# p_ar <- ggplot(plot_df, aes(x = level, y = AR, color = pathogen_label)) +
#   geom_hline(yintercept = 0, linetype = 2, color = "grey65", linewidth = 0.45) +
#   geom_point(position = position_dodge(width = 0.5), size = 1.45) +
#   geom_errorbar(
#     aes(ymin = AR_lo, ymax = AR_hi),
#     position = position_dodge(width = 0.5),
#     width = 0.25,
#     linewidth = 0.5
#   ) +
#   facet_grid(~ subgroup_label, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = my_cols, drop = FALSE) +
#   labs(x = NULL, y = "Attributable mortality (%)", color = NULL) +
#   theme_my +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "none",
#     panel.spacing.x = unit(0.12, "cm")
#   )
# 
# p_ar <- patchwork::wrap_plots(
#   make_panel_title("(b) Attributable mortality"),
#   p_ar,
#   ncol = 1,
#   heights = c(0.10, 1)
# )
# 
# # =========================================================
# # PART C
# # =========================================================
# 
# p_paf <- ggplot(plot_df, aes(x = level, y = PAF, color = pathogen_label)) +
#   geom_hline(yintercept = 0, linetype = 2, color = "grey65", linewidth = 0.45) +
#   geom_point(position = position_dodge(width = 0.5), size = 1.45) +
#   geom_errorbar(
#     aes(ymin = PAF_lo, ymax = PAF_hi),
#     position = position_dodge(width = 0.5),
#     width = 0.25,
#     linewidth = 0.5
#   ) +
#   facet_grid(~ subgroup_label, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = my_cols, drop = FALSE) +
#   guides(
#     color = guide_legend(
#       nrow = 1,
#       byrow = TRUE,
#       override.aes = list(linewidth = 1.1, size = 3)
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Population attributable\nfraction (%)",
#     color = NULL
#   ) +
#   theme_my +
#   theme(
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.box = "horizontal",
#     legend.text = element_text(family = "Times", size = 12, color = "black"),
#     legend.title = element_blank(),
#     legend.key.width = unit(1.4, "lines"),
#     legend.key.height = unit(0.9, "lines"),
#     legend.spacing.x = unit(0.3, "cm"),
#     strip.text = element_blank(),
#     strip.background = element_blank(),
#     axis.text.x = element_text(
#       family = "Times",
#       size = 12,
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5,
#       color = "black"
#     ),
#     panel.spacing.x = unit(0.12, "cm")
#   )
# 
# # =========================================================
# # COMBINE
# # =========================================================
# 
# final_plot <- patchwork::wrap_plots(
#   p_mort,
#   p_ar,
#   p_paf,
#   ncol = 1,
#   heights = c(1.7, 1, 1)
# )
# 
# final_plot
# 
# ggsave(
#   filename = "output/finial_1.pdf",
#   plot = final_plot,
#   width = 13.5,
#   height = 11,
#   device = cairo_pdf
# )
# 
# 
# 
# 
# 
# 
# 
# # =========================================================
# # FINAL COMBINED FIGURE
# # A: Crude mortality as circles
# # B: Attributable mortality
# # C: Population attributable fraction
# # =========================================================
# 
# rm(list = ls())
# 
# suppressPackageStartupMessages({
#   require(pacman)
#   pacman::p_load(
#     dplyr, tidyr, stringr, ggplot2, patchwork,
#     Cairo, grid, survival, scales
#   )
# })
# 
# df_all <- readRDS("data/att_first28.RData")
# final_df <- readRDS("data/att28_paf_marginal.RData")
# 
# # =========================================================
# # COLOURS
# # =========================================================
# 
# my_cols <- c(
#   "AMR"   = "#b22222",
#   "MDR"   = "#Fcae1e",
#   "CRA"   = "#4f7942",
#   "3GCRE" = "#5f8a8b",
#   "CRE"   = "#9e7bb5"
# )
# 
# my_fill <- c(
#   "AMR"     = unname(my_cols["AMR"]),
#   "Non-AMR" = unname(my_cols["AMR"]),
#   "MDR"     = unname(my_cols["MDR"]),
#   "Non-MDR" = unname(my_cols["MDR"]),
#   "CRA"     = unname(my_cols["CRA"]),
#   "CSA"     = unname(my_cols["CRA"]),
#   "3GCRE"   = unname(my_cols["3GCRE"]),
#   "3GCSE"   = unname(my_cols["3GCRE"]),
#   "CRE"     = unname(my_cols["CRE"]),
#   "CSE"     = unname(my_cols["CRE"])
# )
# 
# # =========================================================
# # THEME
# # =========================================================
# 
# theme_my <- theme_minimal(base_size = 12, base_family = "Times") +
#   theme(
#     text = element_text(family = "Times", size = 12, color = "black"),
#     plot.title = element_text(family = "Times", size = 14, face = "bold", color = "black"),
#     axis.title = element_text(family = "Times", size = 12, color = "black"),
#     axis.text = element_text(family = "Times", size = 12, color = "black"),
#     axis.ticks = element_line(color = "black", linewidth = 0.3),
#     axis.ticks.length = unit(-0.1, "cm"),
#     strip.text = element_text(family = "Times", size = 12, color = "black"),
#     strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.35),
#     panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
#     legend.background = element_rect(fill = "white", color = NA),
#     legend.key = element_rect(fill = "white", color = NA)
#   )
# 
# # =========================================================
# # PART A DATA
# # =========================================================
# 
# get_pair_by_subgroup <- function(df, ris_name, labels_pair, pair_name, subgroup_var) {
#   df %>%
#     mutate(
#       group = case_when(
#         grepl("susceptible$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[1],
#         grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[2],
#         TRUE ~ NA_character_
#       ),
#       level = as.character(.data[[subgroup_var]]),
#       subgroup = subgroup_var
#     ) %>%
#     filter(!is.na(group), !is.na(level), level != "") %>%
#     group_by(subgroup, level, group) %>%
#     summarise(
#       total = n(),
#       deaths = sum(event == 1, na.rm = TRUE),
#       mortality = 100 * deaths / total,
#       .groups = "drop"
#     ) %>%
#     mutate(pair = pair_name)
# }
# 
# get_pair_all_subgroups_A <- function(df, ris_name, labels_pair, pair_name) {
#   bind_rows(
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "infection_types"),
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "country_income"),
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "age_group_new")
#   )
# }
# 
# mort_plot_df <- bind_rows(
#   get_pair_all_subgroups_A(df_all[[9]], "amr", c("Non-AMR", "AMR"), "AMR"),
#   get_pair_all_subgroups_A(df_all[[8]], "mdr", c("Non-MDR", "MDR"), "MDR"),
#   get_pair_all_subgroups_A(df_all[[1]], "aci_car", c("CSA", "CRA"), "CRA"),
#   get_pair_all_subgroups_A(df_all[[2]], "ent_thir", c("3GCSE", "3GCRE"), "3GCRE"),
#   get_pair_all_subgroups_A(df_all[[3]], "ent_car", c("CSE", "CRE"), "CRE")
# ) %>%
#   mutate(
#     subgroup_label = case_when(
#       subgroup == "infection_types" ~ "Infection syndromes",
#       subgroup == "country_income" ~ "World Bank income status",
#       subgroup == "age_group_new" ~ "Age groups"
#     ),
#     level = case_when(
#       level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
#       level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
#       level == "Upper middle income" ~ "Upper middle\nincome",
#       level == "Lower middle income" ~ "Lower middle\nincome",
#       TRUE ~ level
#     ),
#     pair = factor(pair, levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")),
#     group = factor(
#       group,
#       levels = c(
#         "AMR", "Non-AMR",
#         "MDR", "Non-MDR",
#         "CRA", "CSA",
#         "3GCRE", "3GCSE",
#         "CRE", "CSE"
#       )
#     ),
#     alpha_group = case_when(
#       group %in% c("Non-AMR", "Non-MDR", "CSA", "3GCSE", "CSE") ~ 0.35,
#       TRUE ~ 0.90
#     )
#   )
# 
# mort_plot_df$level <- factor(
#   mort_plot_df$level,
#   levels = c(
#     "VAP",
#     "Hospital-\nacquired BSI",
#     "Healthcare-\nassociated BSI",
#     "High income",
#     "Upper middle\nincome",
#     "Lower middle\nincome",
#     "<1 year",
#     "1–4 years",
#     "5–14 years",
#     "15–49 years",
#     "50–69 years",
#     "≥70 years"
#   )
# )
# 
# # =========================================================
# # PART A PLOT
# # =========================================================
# 
# make_block_title <- function(title) {
#   ggplot() +
#     annotate(
#       "rect",
#       xmin = -Inf, xmax = Inf,
#       ymin = -Inf, ymax = Inf,
#       fill = "grey90",
#       color = "black",
#       linewidth = 0.35
#     ) +
#     annotate(
#       "text",
#       x = 0.5, y = 0.5,
#       label = title,
#       family = "Times",
#       size = 4
#     ) +
#     xlim(0, 1) +
#     ylim(0, 1) +
#     theme_void() +
#     theme(plot.margin = margin(t = 0, r = 0, b = -5, l = 0))
# }
# 
# make_panel_title <- function(title) {
#   ggplot() +
#     annotate(
#       "text",
#       x = -0.12, y = 0.5,
#       label = title,
#       family = "Times",
#       fontface = "bold",
#       size = 4.2,
#       hjust = 0
#     ) +
#     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
#     theme_void() +
#     theme(plot.margin = margin(t = 6, r = 0, b = 1, l = 0))
# }
# 
# plot_A_block <- function(dat, show_y = FALSE) {
# 
#   p <- ggplot(
#     dat,
#     aes(
#       x = pair,
#       y = mortality,
#       fill = group,
#       alpha = alpha_group
#     )
#   ) +
#     geom_point(
#       position = position_dodge(width = 0.55),
#       shape = 21,
#       size = 3.3,
#       stroke = 0.35,
#       color = "grey35"
#     ) +
#     scale_alpha_identity() +
#     facet_grid(~ level, scales = "free_x", space = "free_x") +
#     scale_fill_manual(values = my_fill, drop = FALSE) +
#     scale_y_continuous(
#       limits = c(0, 90),
#       breaks = seq(0, 80, 20),
#       labels = function(x) paste0(x, "%"),
#       expand = expansion(mult = c(0, 0.04))
#     ) +
#     labs(x = NULL, fill = NULL) +
#     theme_my +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       strip.text.x = element_text(family = "Times", size = 12, color = "black"),
#       panel.spacing.x = unit(0.027, "cm"),
#       plot.margin = margin(t = 0, r = 2, b = 0, l = 2)
#     )
# 
#   if (show_y) {
#     p <- p +
#       labs(y = "Crude 28-day\nmortality (%)") +
#       theme(
#         axis.title.y = element_text(family = "Times", size = 12, color = "black"),
#         axis.text.y = element_text(family = "Times", size = 12, color = "black"),
#         axis.ticks.y = element_line(color = "black", linewidth = 0.25)
#       )
#   } else {
#     p <- p +
#       labs(y = NULL) +
#       theme(
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#       )
#   }
# 
#   p
# }
# 
# p_A_syndrome_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "Infection syndromes"),
#   show_y = TRUE
# )
# 
# p_A_income_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "World Bank income status"),
#   show_y = FALSE
# )
# 
# p_A_age_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "Age groups"),
#   show_y = TRUE
# )
# 
# p_A_syndrome <- patchwork::wrap_plots(
#   make_block_title("Infection syndromes"),
#   p_A_syndrome_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_A_income <- patchwork::wrap_plots(
#   make_block_title("World Bank income status"),
#   p_A_income_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_A_age <- patchwork::wrap_plots(
#   make_block_title("Age groups"),
#   p_A_age_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_mort_body <- (p_A_syndrome | p_A_income) /
#   patchwork::plot_spacer() /
#   p_A_age +
#   plot_layout(heights = c(1, 0.002, 1.05))
# 
# p_mort <- patchwork::wrap_plots(
#   make_panel_title("(a) Crude mortality by resistant–susceptible pairs"),
#   p_mort_body,
#   ncol = 1,
#   heights = c(0.08, 1)
# )
# 
# # =========================================================
# # PART B/C DATA
# # =========================================================
# 
# parse_ci <- function(x) {
#   m <- stringr::str_match(
#     x,
#     "^\\s*([-]?[0-9.]+)%\\s*\\[\\s*([-]?[0-9.]+)\\s*,\\s*([-]?[0-9.]+)\\s*\\]\\s*$"
#   )
#   data.frame(
#     est = as.numeric(m[, 2]),
#     lo  = as.numeric(m[, 3]),
#     hi  = as.numeric(m[, 4])
#   )
# }
# 
# if (!all(c("AR", "AR_lo", "AR_hi") %in% names(final_df))) {
#   tmp <- parse_ci(final_df$AR_CI)
#   final_df$AR <- tmp$est
#   final_df$AR_lo <- tmp$lo
#   final_df$AR_hi <- tmp$hi
# }
# 
# if (!all(c("PAF", "PAF_lo", "PAF_hi") %in% names(final_df))) {
#   tmp <- parse_ci(final_df$PAF_CI)
#   final_df$PAF <- tmp$est
#   final_df$PAF_lo <- tmp$lo
#   final_df$PAF_hi <- tmp$hi
# }
# 
# pathogen_map <- c(
#   aci_car = "CRA",
#   ent_thir = "3GCRE",
#   ent_car = "CRE",
#   pse_car = "CRP",
#   entc_van = "VRE",
#   sa_meth = "MRSA",
#   mdr = "MDR",
#   amr = "AMR"
# )
# 
# plot_df <- final_df %>%
#   mutate(pathogen_label = recode(pathogen, !!!pathogen_map)) %>%
#   filter(
#     pathogen_label %in% c("AMR", "MDR", "CRA", "3GCRE", "CRE"),
#     subgroup %in% c("infection_types", "age_group_new", "country_income"),
#     !is.na(level)
#   ) %>%
#   mutate(
#     subgroup_label = case_when(
#       subgroup == "infection_types" ~ "Infection syndromes",
#       subgroup == "age_group_new" ~ "Age groups",
#       subgroup == "country_income" ~ "World Bank income status"
#     ),
#     level = case_when(
#       level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
#       level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
#       level == "Upper middle income" ~ "Upper middle\nincome",
#       level == "Lower middle income" ~ "Lower middle\nincome",
#       TRUE ~ level
#     ),
#     level = factor(
#       level,
#       levels = c(
#         "VAP",
#         "Hospital-\nacquired BSI",
#         "Healthcare-\nassociated BSI",
#         "<1 year",
#         "1–4 years",
#         "5–14 years",
#         "15–49 years",
#         "50–69 years",
#         "≥70 years",
#         "High income",
#         "Upper middle\nincome",
#         "Lower middle\nincome"
#       )
#     ),
#     subgroup_label = factor(
#       subgroup_label,
#       levels = c(
#         "Infection syndromes",
#         "Age groups",
#         "World Bank income status"
#       )
#     ),
#     pathogen_label = factor(
#       pathogen_label,
#       levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")
#     )
#   )
# 
# # =========================================================
# # PART B
# # =========================================================
# 
# p_ar <- ggplot(plot_df, aes(x = level, y = AR, color = pathogen_label)) +
#   geom_hline(
#     yintercept = 0,
#     linetype = 2,
#     color = "grey65",
#     linewidth = 0.45
#   ) +
#   geom_point(
#     position = position_dodge(width = 0.5),
#     size = 1.45
#   ) +
#   geom_errorbar(
#     aes(ymin = AR_lo, ymax = AR_hi),
#     position = position_dodge(width = 0.5),
#     width = 0.25,
#     linewidth = 0.5
#   ) +
#   facet_grid(~ subgroup_label, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = my_cols, drop = FALSE) +
#   labs(x = NULL, y = "Attributable mortality (%)", color = NULL) +
#   theme_my +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "none",
#     panel.spacing.x = unit(0.12, "cm")
#   )
# 
# p_ar <- patchwork::wrap_plots(
#   make_panel_title("(b) Attributable mortality"),
#   p_ar,
#   ncol = 1,
#   heights = c(0.10, 1)
# )
# 
# # =========================================================
# # PART C
# # =========================================================
# 
# p_paf <- ggplot(plot_df, aes(x = level, y = PAF, color = pathogen_label)) +
#   geom_hline(
#     yintercept = 0,
#     linetype = 2,
#     color = "grey65",
#     linewidth = 0.45
#   ) +
#   geom_point(
#     position = position_dodge(width = 0.5),
#     size = 1.45
#   ) +
#   geom_errorbar(
#     aes(ymin = PAF_lo, ymax = PAF_hi),
#     position = position_dodge(width = 0.5),
#     width = 0.25,
#     linewidth = 0.5
#   ) +
#   facet_grid(~ subgroup_label, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = my_cols, drop = FALSE) +
#   guides(
#     color = guide_legend(
#       nrow = 1,
#       byrow = TRUE,
#       override.aes = list(linewidth = 1.1, size = 3)
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Population attributable\nfraction (%)",
#     color = NULL
#   ) +
#   theme_my +
#   theme(
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.box = "horizontal",
#     legend.text = element_text(family = "Times", size = 12, color = "black"),
#     legend.title = element_blank(),
#     legend.key.width = unit(1.4, "lines"),
#     legend.key.height = unit(0.9, "lines"),
#     legend.spacing.x = unit(0.3, "cm"),
#     strip.text = element_blank(),
#     strip.background = element_blank(),
#     axis.text.x = element_text(
#       family = "Times",
#       size = 12,
#       angle = 90,
#       hjust = 1,
#       vjust = 0.5,
#       color = "black"
#     ),
#     panel.spacing.x = unit(0.12, "cm")
#   )
# 
# # =========================================================
# # COMBINE
# # =========================================================
# 
# final_plot <- patchwork::wrap_plots(
#   p_mort,
#   p_ar,
#   p_paf,
#   ncol = 1,
#   heights = c(1.7, 1, 1)
# )
# 
# final_plot
# 
# ggsave(
#   filename = "output/final_2.pdf",
#   plot = final_plot,
#   width = 13.5,
#   height = 11,
#   device = cairo_pdf
# )


# 
# 
# 
# 
# 
# 
# 
# #####
#Bar
# 
# # =========================================================
# # FINAL COMBINED FIGURE
# # A: Crude mortality
# # B: Attributable mortality
# # C: Population attributable fraction
# # A: separate panel title + block titles
# # =========================================================
# 
# rm(list = ls())
# 
# suppressPackageStartupMessages({
#   require(pacman)
#   pacman::p_load(
#     dplyr, tidyr, stringr, ggplot2, patchwork,
#     Cairo, grid, survival, scales
#   )
# })
# 
# df_all <- readRDS("data/att_first28.RData")
# final_df <- readRDS("data/att28_paf_marginal.RData")
# 
# # =========================================================
# # COLOURS
# # =========================================================
# 
# my_cols <- c(
#   "AMR"   = "#b22222",
#   "MDR"   = "#Fcae1e",
#   "CRA"   = "#4f7942",
#   "3GCRE" = "#5f8a8b",
#   "CRE"   = "#9e7bb5"
# )
# 
# my_fill <- c(
#   "AMR"     = unname(my_cols["AMR"]),
#   "Non-AMR" = scales::alpha(unname(my_cols["AMR"]), 0.35),
#   "MDR"     = unname(my_cols["MDR"]),
#   "Non-MDR" = scales::alpha(unname(my_cols["MDR"]), 0.35),
#   "CRA"     = unname(my_cols["CRA"]),
#   "CSA"     = scales::alpha(unname(my_cols["CRA"]), 0.35),
#   "3GCRE"   = unname(my_cols["3GCRE"]),
#   "3GCSE"   = scales::alpha(unname(my_cols["3GCRE"]), 0.35),
#   "CRE"     = unname(my_cols["CRE"]),
#   "CSE"     = scales::alpha(unname(my_cols["CRE"]), 0.35)
# )
# 
# # =========================================================
# # THEME
# # =========================================================
# 
# theme_my <- theme_minimal(base_size = 12, base_family = "Times") +
#   theme(
#     text = element_text(family = "Times", size = 12, color = "black"),
#     plot.title = element_text(family = "Times", size = 14, face = "bold", color = "black"),
#     axis.title = element_text(family = "Times", size = 12, color = "black"),
#     axis.text = element_text(family = "Times", size = 12, color = "black"),
#     axis.ticks = element_line(color = "black", linewidth = 0.3),
#     axis.ticks.length = unit(-0.1, "cm"),
#     strip.text = element_text(family = "Times", size = 12, color = "black"),
#     strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.35),
#     panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = "white", color = NA),
#     plot.background = element_rect(fill = "white", color = NA),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
#     legend.background = element_rect(fill = "white", color = NA),
#     legend.key = element_rect(fill = "white", color = NA)
#   )
# 
# # =========================================================
# # PART A DATA
# # =========================================================
# 
# get_pair_by_subgroup <- function(df, ris_name, labels_pair, pair_name, subgroup_var) {
#   df %>%
#     mutate(
#       group = case_when(
#         grepl("susceptible$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[1],
#         grepl("resistant$", .data[[ris_name]], ignore.case = TRUE) ~ labels_pair[2],
#         TRUE ~ NA_character_
#       ),
#       level = as.character(.data[[subgroup_var]]),
#       subgroup = subgroup_var
#     ) %>%
#     filter(!is.na(group), !is.na(level), level != "") %>%
#     group_by(subgroup, level, group) %>%
#     summarise(
#       total = n(),
#       deaths = sum(event == 1, na.rm = TRUE),
#       mortality = 100 * deaths / total,
#       .groups = "drop"
#     ) %>%
#     mutate(pair = pair_name)
# }
# 
# get_pair_all_subgroups_A <- function(df, ris_name, labels_pair, pair_name) {
#   bind_rows(
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "infection_types"),
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "country_income"),
#     get_pair_by_subgroup(df, ris_name, labels_pair, pair_name, "age_group_new")
#   )
# }
# 
# mort_plot_df <- bind_rows(
#   get_pair_all_subgroups_A(df_all[[9]], "amr", c("Non-AMR", "AMR"), "AMR"),
#   get_pair_all_subgroups_A(df_all[[8]], "mdr", c("Non-MDR", "MDR"), "MDR"),
#   get_pair_all_subgroups_A(df_all[[1]], "aci_car", c("CSA", "CRA"), "CRA"),
#   get_pair_all_subgroups_A(df_all[[2]], "ent_thir", c("3GCSE", "3GCRE"), "3GCRE"),
#   get_pair_all_subgroups_A(df_all[[3]], "ent_car", c("CSE", "CRE"), "CRE")
# ) %>%
#   mutate(
#     subgroup_label = case_when(
#       subgroup == "infection_types" ~ "Infection syndromes",
#       subgroup == "country_income" ~ "World Bank income status",
#       subgroup == "age_group_new" ~ "Age groups"
#     ),
#     level = case_when(
#       level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
#       level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
#       level == "Upper middle income" ~ "Upper middle\nincome",
#       level == "Lower middle income" ~ "Lower middle\nincome",
#       TRUE ~ level
#     ),
#     pair = factor(pair, levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")),
#     group = factor(
#       group,
#       levels = c(
#         "AMR", "Non-AMR",
#         "MDR", "Non-MDR",
#         "CRA", "CSA",
#         "3GCRE", "3GCSE",
#         "CRE", "CSE"
#       )
#     )
#   )
# 
# mort_plot_df$level <- factor(
#   mort_plot_df$level,
#   levels = c(
#     "VAP",
#     "Hospital-\nacquired BSI",
#     "Healthcare-\nassociated BSI",
#     "High income",
#     "Upper middle\nincome",
#     "Lower middle\nincome",
#     "<1 year",
#     "1–4 years",
#     "5–14 years",
#     "15–49 years",
#     "50–69 years",
#     "≥70 years"
#   )
# )
# 
# # =========================================================
# # PART A PLOT
# # =========================================================
# 
# make_block_title <- function(title) {
#   ggplot() +
#     annotate(
#       "rect",
#       xmin = -Inf,
#       xmax = Inf,
#       ymin = -Inf,
#       ymax = Inf,
#       fill = "grey90",
#       color = "black",
#       linewidth = 0.35
#     ) +
#     annotate(
#       "text",
#       x = 0.5,
#       y = 0.5,
#       label = title,
#       family = "Times",
#       size = 4
#     ) +
#     xlim(0, 1) +
#     ylim(0, 1) +
#     theme_void() +
#     theme(
#       plot.margin = margin(t = 0, r = 0, b = -5, l = 0)
#     )
# }
# 
# make_panel_title <- function(title) {
#   ggplot() +
#     annotate(
#       "text",
#       x = -0.12,
#       y = 0.5,
#       label = title,
#       family = "Times",
#       fontface = "bold",
#       size = 4.2,
#       hjust = 0
#     ) +
#     coord_cartesian(
#       xlim = c(0, 1),
#       ylim = c(0, 1),
#       clip = "off"
#     ) +
#     theme_void() +
#     theme(
#       plot.margin = margin(t = 6, r = 0, b = 1, l = 0)
#     )
# }
# 
# plot_A_block <- function(dat, show_y = FALSE) {
#   
#   p <- ggplot(dat, aes(x = pair, y = mortality, fill = group)) +
#     geom_col(
#       position = position_dodge(width = 0.7),
#       width = 0.55,
#       color = "black",
#       linewidth = 0.3
#     ) +
#     facet_grid(~ level, scales = "free_x", space = "free_x") +
#     scale_fill_manual(values = my_fill, drop = FALSE) +
#     scale_y_continuous(
#       limits = c(0, 90),
#       breaks = seq(0, 80, 20),
#       expand = expansion(mult = c(0, 0.04))
#     ) +
#     labs(x = NULL, fill = NULL) +
#     theme_my +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       strip.text.x = element_text(
#         family = "Times", size = 12, color = "black"
#       ),
#       panel.spacing.x = unit(0.027, "cm"),
#       plot.margin = margin(t = 0, r = 2, b = 0, l = 2)
#     )
#   
#   if (show_y) {
#     p <- p +
#       labs(y = "Crude 28-day\nmortality (%)") +
#       theme(
#         axis.title.y = element_text(family = "Times", size = 12, color = "black"),
#         axis.text.y = element_text(family = "Times", size = 12, color = "black"),
#         axis.ticks.y = element_line(color = "black", linewidth = 0.25)
#       )
#   } else {
#     p <- p +
#       labs(y = NULL) +
#       theme(
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()
#       )
#   }
#   
#   p
# }
# 
# p_A_syndrome_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "Infection syndromes"),
#   show_y = TRUE
# )
# 
# p_A_income_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "World Bank income status"),
#   show_y = FALSE
# )
# 
# p_A_age_body <- plot_A_block(
#   mort_plot_df %>% filter(subgroup_label == "Age groups"),
#   show_y = TRUE
# )
# 
# p_A_syndrome <- patchwork::wrap_plots(
#   make_block_title("Infection syndromes"),
#   p_A_syndrome_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_A_income <- patchwork::wrap_plots(
#   make_block_title("World Bank income status"),
#   p_A_income_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# p_A_age <- patchwork::wrap_plots(
#   make_block_title("Age groups"),
#   p_A_age_body,
#   ncol = 1,
#   heights = c(0.3, 1)
# )
# 
# layout_A <- c(
#   area(t = 1, l = 1, b = 1, r = 1),  # Infection
#   area(t = 1, l = 2, b = 1, r = 2),  # Income
#   area(t = 2, l = 1, b = 2, r = 2)   # Age spans two columns
# )
# 
# p_mort_body <- (p_A_syndrome | p_A_income) /
#   patchwork::plot_spacer() /
#   p_A_age +
#   plot_layout(
#     heights = c(1, 0.002, 1.05)
#   )
# 
# p_mort <- patchwork::wrap_plots(
#   make_panel_title("(a) Crude mortality by resistant–susceptible pairs"),
#   p_mort_body,
#   ncol = 1,
#   heights = c(0.08, 1)
# )
# 
# # =========================================================
# # PART B/C DATA
# # =========================================================
# 
# parse_ci <- function(x) {
#   m <- stringr::str_match(
#     x,
#     "^\\s*([-]?[0-9.]+)%\\s*\\[\\s*([-]?[0-9.]+)\\s*,\\s*([-]?[0-9.]+)\\s*\\]\\s*$"
#   )
#   data.frame(
#     est = as.numeric(m[, 2]),
#     lo  = as.numeric(m[, 3]),
#     hi  = as.numeric(m[, 4])
#   )
# }
# 
# if (!all(c("AR", "AR_lo", "AR_hi") %in% names(final_df))) {
#   tmp <- parse_ci(final_df$AR_CI)
#   final_df$AR <- tmp$est
#   final_df$AR_lo <- tmp$lo
#   final_df$AR_hi <- tmp$hi
# }
# 
# if (!all(c("PAF", "PAF_lo", "PAF_hi") %in% names(final_df))) {
#   tmp <- parse_ci(final_df$PAF_CI)
#   final_df$PAF <- tmp$est
#   final_df$PAF_lo <- tmp$lo
#   final_df$PAF_hi <- tmp$hi
# }
# 
# pathogen_map <- c(
#   aci_car = "CRA",
#   ent_thir = "3GCRE",
#   ent_car = "CRE",
#   pse_car = "CRP",
#   entc_van = "VRE",
#   sa_meth = "MRSA",
#   mdr = "MDR",
#   amr = "AMR"
# )
# 
# plot_df <- final_df %>%
#   mutate(pathogen_label = recode(pathogen, !!!pathogen_map)) %>%
#   filter(
#     pathogen_label %in% c("AMR", "MDR", "CRA", "3GCRE", "CRE"),
#     subgroup %in% c("infection_types", "age_group_new", "country_income"),
#     !is.na(level)
#   ) %>%
#   mutate(
#     subgroup_label = case_when(
#       subgroup == "infection_types" ~ "Infection syndromes",
#       subgroup == "age_group_new" ~ "Age groups",
#       subgroup == "country_income" ~ "World Bank income status"
#     ),
#     level = case_when(
#       level == "Hospital-acquired BSI" ~ "Hospital-\nacquired BSI",
#       level == "Healthcare-associated BSI" ~ "Healthcare-\nassociated BSI",
#       level == "Upper middle income" ~ "Upper middle\nincome",
#       level == "Lower middle income" ~ "Lower middle\nincome",
#       TRUE ~ level
#     ),
#     level = factor(
#       level,
#       levels = c(
#         "VAP",
#         "Hospital-\nacquired BSI",
#         "Healthcare-\nassociated BSI",
#         "<1 year",
#         "1–4 years",
#         "5–14 years",
#         "15–49 years",
#         "50–69 years",
#         "≥70 years",
#         "High income",
#         "Upper middle\nincome",
#         "Lower middle\nincome"
#       )
#     ),
#     subgroup_label = factor(
#       subgroup_label,
#       levels = c(
#         "Infection syndromes",
#         "Age groups",
#         "World Bank income status"
#       )
#     ),
#     pathogen_label = factor(
#       pathogen_label,
#       levels = c("AMR", "MDR", "CRA", "3GCRE", "CRE")
#     )
#   )
# 
# # =========================================================
# # PART B
# # =========================================================
# 
# p_ar <- ggplot(plot_df, aes(x = level, y = AR, color = pathogen_label)) +
#   geom_hline(
#     yintercept = 0,
#     linetype = 2,
#     color = "grey65",
#     linewidth = 0.45
#   ) +
#   geom_point(
#     position = position_dodge(width = 0.5),
#     size = 1.45
#   ) +
#   geom_errorbar(
#     aes(ymin = AR_lo, ymax = AR_hi),
#     position = position_dodge(width = 0.5),
#     width = 0.25,
#     linewidth = 0.5
#   ) +
#   facet_grid(~ subgroup_label, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = my_cols, drop = FALSE) +
#   labs(x = NULL, y = "Attributable mortality (%)", color = NULL) +
#   theme_my +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "none",
#     panel.spacing.x = unit(0.12, "cm")
#   ) 
# 
# p_ar_body <- p_ar
# 
# p_ar <- patchwork::wrap_plots(
#   make_panel_title("(b) Attributable mortality"),
#   p_ar_body,
#   ncol = 1,
#   heights = c(0.10, 1)
# )
# 
# # =========================================================
# # PART C
# # =========================================================
# 
# p_paf <- ggplot(plot_df, aes(x = level, y = PAF, color = pathogen_label)) +
#   geom_hline(
#     yintercept = 0,
#     linetype = 2,
#     color = "grey65",
#     linewidth = 0.45
#   ) +
#   geom_point(
#     position = position_dodge(width = 0.5),
#     size = 1.45
#   ) +
#   geom_errorbar(
#     aes(ymin = PAF_lo, ymax = PAF_hi),
#     position = position_dodge(width = 0.5),
#     width = 0.25,
#     linewidth = 0.5
#   ) +
#   facet_grid(~ subgroup_label, scales = "free_x", space = "free_x") +
#   scale_color_manual(values = my_cols, drop = FALSE) +
#   guides(
#     color = guide_legend(
#       nrow = 1,
#       byrow = TRUE,
#       override.aes = list(
#         linewidth = 1.1,
#         size = 3
#       )
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Population attributable\nfraction (%)",
#     color = NULL
#   ) +
#   theme_my +
#   theme(
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.box = "horizontal",
#     legend.text = element_text(family = "Times", size = 12, color = "black"),
#     legend.title = element_blank(),
#     legend.key.width = unit(1.4, "lines"),
#     legend.key.height = unit(0.9, "lines"),
#     legend.spacing.x = unit(0.3, "cm"),
#     strip.text = element_blank(),
#     strip.background = element_blank(),
#     axis.text.x = element_text(
#       family = "Times", size = 12, angle = 90,
#       hjust = 1, vjust = 0.5, color = "black"
#     ),
#     panel.spacing.x = unit(0.12, "cm")
#   )
# 
# # =========================================================
# # COMBINE
# # =========================================================
# 
# final_plot <- patchwork::wrap_plots(
#   p_mort,
#   p_ar,
#   p_paf,
#   ncol = 1,
#   heights = c(1.7, 1, 1)
# )
# 
# final_plot
# 
# ggsave(
#   filename = "output/final_combined_figure_revised.pdf",
#   plot = final_plot,
#   width = 12,
#   height = 11,
#   device = cairo_pdf
# )
# 
