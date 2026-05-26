# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    ggplot2,
    scales,
    patchwork,
    Cairo,
    extrafont,
    viridisLite,
    grid,
    ggtext
  )
})

# Working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df_plot_all  <- readRDS("data/data_plot_all.RData")
df_count_all <- readRDS("data/df_count_all.RData")

# Basic cleaning
df_plot_all <- df_plot_all %>%
  mutate(
    pathogen   = as.character(pathogen),
    subgroup   = as.character(subgroup),
    level      = as.character(level),
    group      = as.character(group),
    fbis_score = factor(as.character(fbis_score), levels = as.character(1:7), ordered = TRUE)
  )

df_count_all <- df_count_all %>%
  mutate(
    pathogen   = as.character(pathogen),
    subgroup   = as.character(subgroup),
    level      = as.character(level),
    group      = as.character(group),
    fbis_score = factor(as.character(fbis_score), levels = as.character(1:7), ordered = TRUE)
  )

# Standardise resistant / susceptible
standardise_group <- function(x) {
  case_when(
    grepl("res", x, ignore.case = TRUE) ~ "resistant",
    grepl("sus", x, ignore.case = TRUE) ~ "susceptible",
    TRUE ~ NA_character_
  )
}

df_plot_all <- df_plot_all %>%
  mutate(group_std = standardise_group(group))

df_count_all <- df_count_all %>%
  mutate(group_std = standardise_group(group))

# ============================================================
# FBIS labels / heatmap
# ============================================================
fbis_labels <- data.frame(
  fbis_score = factor(1:7, levels = as.character(1:7), ordered = TRUE),
  severity = 1:7,
  label = c(
    "On palliative care in terminal phases.",
    "Accommodated in a long-term ventilator unit.",
    "Hospitalized in an intensive care unit.",
    "Hospitalized but not requiring\nan intensive care unit.",
    "Out of hospital; significant disability;\nrequires assistance.",
    "Out of hospital; moderate signs or symptoms;\nunable to complete daily activities.",
    "Out of hospital; basically healthy;\nable to complete activities."
  )
)

p_heatmap <- ggplot(fbis_labels, aes(x = 1, y = fbis_score, fill = severity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = label),
    color = "white",
    fontface = "bold",
    size = 3.8,
    hjust = 0.5
  ) +
  scale_fill_gradientn(
    colors = viridisLite::viridis(7, option = "plasma", end = 0.9, direction = 1),
    guide = "none"
  ) +
  labs(
    title = "FBIS score",
    x = NULL,
    y = NULL
  ) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(6, "mm"),
    plot.title = element_text(
      hjust = 0.5,
      family = "Times New Roman",
      size = 12,
      face = "bold"
    ),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================================
# Orders / labels
# ============================================================
pathogen_order <- c("aci_car", "ent_thir", "ent_car", "pse_car", "entc_van", "sa_meth")

pathogen_titles <- c(
  aci_car  = "CRA versus CSA",
  ent_thir = "3GCRE versus 3GCSE",
  ent_car  = "CRE versus CSE",
  pse_car  = "CRP versus CSP",
  entc_van = "VRE versus VSE",
  sa_meth  = "MRSA versus MSSA"
)

pathogen_title_map_md <- c(
  aci_car  = "<i>Acinetobacter</i> spp.",
  ent_thir = "Enterobacterales",
  ent_car  = "Enterobacterales",
  pse_car  = "<i>Pseudomonas aeruginosa</i>",
  entc_van = "<i>Enterococcus</i> spp.",
  sa_meth  = "<i>Staphylococcus aureus</i>"
)

pretty_group_labels <- function(pathogen_name, group_std) {
  if (pathogen_name == "aci_car") {
    return(ifelse(group_std == "resistant", "Carbapenem-resistant", "Carbapenem-susceptible"))
  }
  if (pathogen_name == "ent_thir") {
    return(ifelse(
      group_std == "resistant",
      "Third-generation cephalosporin-resistant",
      "Third-generation cephalosporin-susceptible"
    ))
  }
  if (pathogen_name == "ent_car") {
    return(ifelse(group_std == "resistant", "Carbapenem-resistant", "Carbapenem-susceptible"))
  }
  if (pathogen_name == "pse_car") {
    return(ifelse(group_std == "resistant", "Carbapenem-resistant", "Carbapenem-susceptible"))
  }
  if (pathogen_name == "entc_van") {
    return(ifelse(group_std == "resistant", "Vancomycin-resistant", "Vancomycin-susceptible"))
  }
  if (pathogen_name == "sa_meth") {
    return(ifelse(group_std == "resistant", "Methicillin-resistant", "Methicillin-susceptible"))
  }
  group_std
}

# ============================================================
# Order levels by subgroup
# ============================================================
order_levels_by_subgroup <- function(x, subgroup_name) {
  x <- as.character(x)
  
  if (subgroup_name == "overall") {
    lev <- "overall"
    return(factor(x, levels = lev))
  }
  
  if (subgroup_name == "infection_types") {
    preferred <- c(
      "VAP",
      "Hospital-acquired BSI",
      "Healthcare-associated BSI"
    )
    lev <- c(preferred[preferred %in% x], setdiff(unique(x), preferred))
    return(factor(x, levels = lev))
  }
  
  if (subgroup_name == "country_income") {
    preferred <- c(
      "High income",
      "Upper middle income",
      "Lower middle income",
      "Low income"
    )
    lev <- c(preferred[preferred %in% x], setdiff(unique(x), preferred))
    return(factor(x, levels = lev))
  }
  
  factor(x)
}

# ============================================================
# OVERALL SECTION
# ============================================================

# ------------------------------------------------------------
# Merge plot + count for one pathogen/subgroup/level
# ------------------------------------------------------------
get_panel_data <- function(pathogen_name, subgroup_name, level_name) {
  
  pred <- df_plot_all %>%
    dplyr::filter(
      pathogen == pathogen_name,
      subgroup == subgroup_name,
      level == level_name
    )
  
  cnt <- df_count_all %>%
    dplyr::filter(
      pathogen == pathogen_name,
      subgroup == subgroup_name,
      level == level_name
    ) %>%
    dplyr::select(fbis_score, group_std, expected_count, N)
  
  if (nrow(pred) == 0) return(NULL)
  
  out <- pred %>%
    left_join(cnt, by = c("fbis_score", "group_std")) %>%
    mutate(
      group = pretty_group_labels(pathogen_name, group_std)
    )
  
  out
}

# ------------------------------------------------------------
# Make one OVERALL panel
# ------------------------------------------------------------
make_one_panel <- function(pathogen_name, subgroup_name, level_name,
                           show_x = FALSE, show_y = FALSE, show_legend = FALSE) {
  
  dat <- get_panel_data(pathogen_name, subgroup_name, level_name)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  
  ggplot(dat, aes(x = factor(fbis_score),
                  y = mean_prob * 100,
                  color = group, group = group)) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 0.5) +
    geom_errorbar(
      aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      width = 0.3, linewidth = 0.5
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      labels = function(x) ifelse(x %% 20 == 0, scales::label_percent(scale = 1)(x), ""),
      limits = c(0, 100)
    ) +
    labs(
      title = pathogen_titles[[pathogen_name]],
      x = if (show_x) "FBIS score" else "",
      y = if (show_y) "Predicted probability" else ""
    ) +
    scale_color_manual(
      name = " ",
      values = c("#0073c2FF", "#EFC000FF"),
      labels = c("Resistant (R)", "Susceptible (S)")
    ) +
    geom_text(
      data = dat %>% dplyr::filter(group_std == "susceptible"),
      aes(x = fbis_score, y = 79, label = paste0("S = ", expected_count)),
      color = "#EFC000FF",
      size = 3.6,
      family = "Times New Roman",
      hjust = 0,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = dat %>% dplyr::filter(group_std == "resistant"),
      aes(x = fbis_score, y = 61, label = paste0("R = ", expected_count)),
      color = "#0073c2FF",
      size = 3.6,
      family = "Times New Roman",
      hjust = 0,
      inherit.aes = FALSE
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Times New Roman", size = 12),
      strip.text = element_text(family = "Times New Roman", size = 12),
      axis.text.y = if (show_y) {
        element_text(size = 12, family = "Times New Roman", color = "black")
      } else {
        element_blank()
      },
      axis.text.x = element_text(
        size = 12, family = "Times New Roman",
        color = "black", margin = margin(t = 3)
      ),
      axis.title.x = if (show_x) {
        element_text(family = "Times New Roman", size = 12, margin = margin(t = 5))
      } else {
        element_blank()
      },
      axis.title.y = if (show_y) {
        element_text(family = "Times New Roman", size = 12, margin = margin(r = 8))
      } else {
        element_blank()
      },
      plot.title = element_text(
        family = "Times New Roman",
        size = 12, face = "bold", hjust = 0.5
      ),
      legend.title = element_text(size = 12, family = "Times New Roman"),
      legend.text = element_text(size = 12, family = "Times New Roman"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.margin = margin(t = -7, r = 0, b = 0, l = 0),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
    ) +
    coord_flip()
}
# ------------------------------------------------------------
# Build OVERALL figure
# ------------------------------------------------------------
make_overall_figure <- function(file_prefix) {
  
  panel_list <- vector("list", 6)
  
  for (i in seq_along(pathogen_order)) {
    panel_list[[i]] <- make_one_panel(
      pathogen_name = pathogen_order[i],
      subgroup_name = "overall",
      level_name    = "overall",
      show_x        = (i %in% c(2, 5)),
      show_y        = (i %in% c(1, 4)),
      show_legend   = (i %in% c(2, 5))
    )
  }
  
  if (all(sapply(panel_list, is.null))) return(NULL)
  
  # Row 1
  final_plot1 <- p_heatmap + panel_list[[1]] + panel_list[[2]] + panel_list[[3]] +
    plot_layout(widths = c(2, 2, 2, 2))
  
  # Row 2
  final_plot2 <- p_heatmap + panel_list[[4]] + panel_list[[5]] + panel_list[[6]] +
    plot_layout(widths = c(2, 2, 2, 2))
  
  # Save row 1
  CairoPDF(paste0(file_prefix, "1.pdf"), width = 16, height = 7)
  print(final_plot1)
  dev.off()
  
  # Save row 2
  CairoPDF(paste0(file_prefix, "2.pdf"), width = 16, height = 7)
  print(final_plot2)
  dev.off()
  
  # ------------------------------------------------------------
  # Combined figure (single legend)
  # ------------------------------------------------------------
  
  row1_noleg <- p_heatmap +
    make_one_panel("aci_car",  "overall", "overall", show_x = FALSE, show_y = TRUE,  show_legend = FALSE) +
    make_one_panel("ent_thir", "overall", "overall", show_x = TRUE,  show_y = FALSE, show_legend = FALSE) +
    make_one_panel("ent_car",  "overall", "overall", show_x = FALSE, show_y = FALSE, show_legend = FALSE) +
    plot_layout(widths = c(2, 2, 2, 2))
  
  row2_withleg <- p_heatmap +
    make_one_panel("pse_car",  "overall", "overall", show_x = FALSE, show_y = TRUE,  show_legend = FALSE) +
    make_one_panel("entc_van", "overall", "overall", show_x = TRUE,  show_y = FALSE, show_legend = TRUE) +
    make_one_panel("sa_meth",  "overall", "overall", show_x = FALSE, show_y = FALSE, show_legend = FALSE) +
    plot_layout(widths = c(2, 2, 2, 2))
  
  final_plot_combined <- (row1_noleg / row2_withleg) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  CairoPDF(paste0(file_prefix, "_combined.pdf"), width = 16, height = 12)
  print(final_plot_combined)
  dev.off()
  
  list(
    plot1 = final_plot1,
    plot2 = final_plot2,
    plot_combined = final_plot_combined
  )
}

# Run OVERALL
make_overall_figure("output/fbis_overall")

# ============================================================
# SUBGROUP SECTION 
# infection_types + country_income only
# age not plotted
# ============================================================

# ------------------------------------------------------------
# Merge plot + count for one pathogen + subgroup
# ------------------------------------------------------------
get_subgroup_plot_data <- function(pathogen_name, subgroup_name) {
  
  pred <- df_plot_all %>%
    dplyr::filter(
      pathogen == pathogen_name,
      subgroup == subgroup_name,
      !is.na(level),
      level != "NA"
    )
  
  cnt <- df_count_all %>%
    dplyr::filter(
      pathogen == pathogen_name,
      subgroup == subgroup_name,
      !is.na(level),
      level != "NA"
    ) %>%
    dplyr::select(pathogen, subgroup, level, fbis_score, group_std, expected_count, N)
  
  if (nrow(pred) == 0) return(NULL)
  
  out <- pred %>%
    left_join(
      cnt,
      by = c("pathogen", "subgroup", "level", "fbis_score", "group_std")
    ) %>%
    mutate(
      level = order_levels_by_subgroup(level, subgroup_name),
      group = case_when(
        group_std == "resistant"   ~ pretty_group_labels(pathogen_name, "resistant"),
        group_std == "susceptible" ~ pretty_group_labels(pathogen_name, "susceptible"),
        TRUE ~ group
      )
    ) %>%
    dplyr::filter(!is.na(level))
  
  out
}

# ------------------------------------------------------------
# Make one subgroup plot
# ------------------------------------------------------------
make_subgroup_plot <- function(pathogen_name, subgroup_name) {
  
  dat <- get_subgroup_plot_data(pathogen_name, subgroup_name)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  
  n_facets <- length(unique(dat$level))
  ncol_use <- ifelse(n_facets <= 3, n_facets, 3)
  
  p <- ggplot(
    dat,
    aes(
      x = factor(fbis_score),
      y = mean_prob * 100,
      color = group_std,
      group = group_std
    )
  ) +
    geom_point(size = 1.5) +
    geom_line(linewidth = 0.5) +
    geom_errorbar(
      aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      width = 0.4,
      linewidth = 0.5
    ) +
    facet_wrap(~ level, ncol = ncol_use) +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      breaks = seq(0, 100, by = 10),
      limits = c(0, 100)
    ) +
    labs(
      title = "",
      x = "FBIS score",
      y = "Predicted probability"
    ) +
    scale_color_manual(
      name = " ",
      values = c(
        resistant   = "#0073c2FF",
        susceptible = "#EFC000FF"
      ),
      breaks = c("resistant", "susceptible"),
      labels = c("Resistant (R)", "Susceptible (S)")
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Times New Roman", size = 10),
      plot.title = element_text(
        family = "Times New Roman",
        size = 10,
        face = "bold",
        hjust = 0.5
      ),
      strip.text = element_text(family = "Times New Roman", size = 10),
      axis.title.x = element_text(
        family = "Times New Roman",
        size = 10,
        margin = margin(t = 5)
      ),
      axis.title.y = element_text(
        family = "Times New Roman",
        size = 10,
        margin = margin(r = 5)
      ),
      axis.text.y = element_text(
        size = 10,
        family = "Times New Roman",
        color = "black",
        margin = margin(r = 10)
      ),
      axis.text.x = element_text(
        size = 10,
        family = "Times New Roman",
        color = "black",
        margin = margin(t = 3)
      ),
      legend.title = element_text(size = 10, family = "Times New Roman"),
      legend.text = element_text(size = 10, family = "Times New Roman"),
      panel.grid.major = element_line(
        color = "gray80",
        linetype = "dotted",
        linewidth = 0.3
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        linewidth = 0.3,
        color = "black",
        fill = NA
      ),
      strip.background = element_rect(
        linewidth = 0.3,
        color = "black",
        fill = "gray80"
      ),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(t = -7, r = 0, b = 0, l = 0),
      plot.margin = margin(t = -10, r = 0, b = 0, l = 0)
    )
  
  p
}

# ------------------------------------------------------------
# Save infection_types plots
# ------------------------------------------------------------
# for (pathogen_name in pathogen_order) {
#   
#   p_inf <- make_subgroup_plot(pathogen_name, "infection_types")
#   
#   if (!is.null(p_inf)) {
#     CairoPDF(paste0("output/", pathogen_name, "_infection_types.pdf"), width = 8, height = 5.5)
#     print(p_inf)
#     dev.off()
#   }
# }

# ------------------------------------------------------------
# Save country_income plots
# ------------------------------------------------------------
# for (pathogen_name in pathogen_order) {
#   
#   p_inc <- make_subgroup_plot(pathogen_name, "country_income")
#   
#   if (!is.null(p_inc)) {
#     CairoPDF(paste0("output/", pathogen_name, "_country_income.pdf"), width = 8, height = 5.5)
#     print(p_inc)
#     dev.off()
#   }
# }


# ============================================================
# COMBINE SUBGROUP PLOTS INTO ONE FIGURE
# infection_types -> one combined figure
# country_income  -> one combined figure
# ============================================================

# ------------------------------------------------------------
# Helper: one compact panel for combined figures
# ------------------------------------------------------------
make_subgroup_panel_compact <- function(pathogen_name, subgroup_name,
                                        show_y = FALSE, show_legend = FALSE) {
  
  dat <- get_subgroup_plot_data(pathogen_name, subgroup_name)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  
  ggplot(
    dat,
    aes(
      x = factor(fbis_score),
      y = mean_prob * 100,
      color = group_std,
      group = group_std
    )
  ) +
    geom_point(size = 1.2) +
    geom_line(linewidth = 0.45) +
    geom_errorbar(
      aes(ymin = lower_ci * 100, ymax = upper_ci * 100),
      width = 0.35,
      linewidth = 0.4
    ) +
    facet_wrap(~ level, ncol = length(unique(dat$level))) +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      breaks = seq(0, 100, by = 20),
      limits = c(0, 100)
    ) +
    labs(
      title = pathogen_titles[[pathogen_name]],
      x = NULL,
      y = if (show_y) "Predicted probability" else ""
    ) +
    scale_color_manual(
      name = " ",
      values = c(
        resistant   = "#0073c2FF",
        susceptible = "#EFC000FF"
      ),
      breaks = c("resistant", "susceptible"),
      labels = c("Resistant (R)", "Susceptible (S)")
    ) +
    theme_bw() +
    theme(
      text = element_text(family = "Times New Roman", size = 9),
      plot.title = element_text(
        family = "Times New Roman",
        size = 10,
        face = "bold",
        hjust = 0.5
      ),
      strip.text = element_text(family = "Times New Roman", size = 9),
      axis.title.x = element_text(
        family = "Times New Roman",
        size = 9,
        margin = margin(t = 4)
      ),
      axis.title.y = element_text(
        family = "Times New Roman",
        size = 9,
        margin = margin(r = 4)
      ),
      axis.text.y = if (show_y) {
        element_text(size = 8.5, family = "Times New Roman", color = "black")
      } else {
        element_blank()
      },
      axis.text.x = element_text(
        size = 8.5,
        family = "Times New Roman",
        color = "black",
        margin = margin(t = 2)
      ),
      legend.title = element_text(size = 9, family = "Times New Roman"),
      legend.text = element_text(size = 8.5, family = "Times New Roman"),
      panel.grid.major = element_line(
        color = "gray85",
        linetype = "dotted",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        linewidth = 0.3,
        color = "black",
        fill = NA
      ),
      strip.background = element_rect(
        linewidth = 0.3,
        color = "black",
        fill = "gray80"
      ),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )
}

# ------------------------------------------------------------
# Build combined figure by subgroup
# 6 pathogens in 2 rows x 3 columns
# ------------------------------------------------------------
make_combined_subgroup_figure <- function(subgroup_name, file_name, plot_title = NULL) {
  
  if (subgroup_name == "infection_types") {
    
    # VRE removed
    p1 <- make_subgroup_panel_compact("aci_car",  subgroup_name, show_y = TRUE,  show_legend = FALSE)
    p2 <- make_subgroup_panel_compact("ent_thir", subgroup_name, show_y = FALSE, show_legend = TRUE)
    p3 <- make_subgroup_panel_compact("ent_car",  subgroup_name, show_y = FALSE, show_legend = FALSE)
    p4 <- make_subgroup_panel_compact("pse_car",  subgroup_name, show_y = TRUE,  show_legend = FALSE)
    p5 <- make_subgroup_panel_compact("sa_meth",  subgroup_name, show_y = FALSE, show_legend = FALSE)
    
    p5 <- p5 + labs(x = "FBIS score")
    
    plot_list <- list(p1, p2, p3, p4, p5)
    
  } else {
    
    # income still all 6
    p1 <- make_subgroup_panel_compact("aci_car",  subgroup_name, show_y = TRUE,  show_legend = FALSE)
    p2 <- make_subgroup_panel_compact("ent_thir", subgroup_name, show_y = FALSE, show_legend = TRUE)
    p3 <- make_subgroup_panel_compact("ent_car",  subgroup_name, show_y = FALSE, show_legend = FALSE)
    p4 <- make_subgroup_panel_compact("pse_car",  subgroup_name, show_y = TRUE,  show_legend = FALSE)
    p5 <- make_subgroup_panel_compact("entc_van", subgroup_name, show_y = FALSE, show_legend = FALSE)
    p6 <- make_subgroup_panel_compact("sa_meth",  subgroup_name, show_y = FALSE, show_legend = FALSE)
    
    p5 <- p5 + labs(x = "FBIS score")
    
    plot_list <- list(p1, p2, p3, p4, p5, p6)
    
  }
  
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  combined_plot <- wrap_plots(plot_list, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  
  if (!is.null(plot_title)) {
    combined_plot <- combined_plot +
      plot_annotation(
        title = plot_title,
        theme = theme(
          plot.title = element_text(
            family = "Times New Roman",
            size = 11,
            face = "bold",
            hjust = 0.5
          )
        )
      )
  }
  
  CairoPDF(file_name, width = 16, height = 9)
  print(combined_plot)
  dev.off()
  
  combined_plot
}

# ------------------------------------------------------------
# Combined infection_types figure
# ------------------------------------------------------------
p_infection_all <- make_combined_subgroup_figure(
  subgroup_name = "infection_types",
  file_name     = "output/fbis_infection_types_combined.pdf",
  plot_title    = ""
)

# ------------------------------------------------------------
# Combined country_income figure
# ------------------------------------------------------------
p_income_all <- make_combined_subgroup_figure(
  subgroup_name = "country_income",
  file_name     = "output/fbis_country_income_combined.pdf",
  plot_title    = ""
)

