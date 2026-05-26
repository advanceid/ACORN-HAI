# =========================================================
# EQ-5D-3L adjusted mean plot
# 2 x 3 panel layout
# one line per subgroup
# R = blue, S = yellow
# =========================================================

rm(list = ls())

suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    ggplot2,
    patchwork,
    Cairo
  )
})

# --------------------------------------------------
# Working directory
# --------------------------------------------------
wd <- "./"
setwd(wd)

# --------------------------------------------------
# Load data
# --------------------------------------------------
df_count_all <- readRDS("data/eq5d3l_count_all.RData") %>%
  as.data.frame()

# --------------------------------------------------
# Labels
# --------------------------------------------------
keep_pathogen <- c(
  "aci_car", "ent_thir", "ent_car",
  "pse_car", "entc_van", "sa_meth"
)

pathogen_map <- c(
  aci_car  = "CRA",
  ent_thir = "3GCRE",
  ent_car  = "CRE",
  pse_car  = "CRP",
  entc_van = "VRE",
  sa_meth  = "MRSA"
)

title_map <- c(
  CRA     = "CRA versus CSA",
  "3GCRE" = "3GCRE versus 3GCSE",
  CRE     = "CRE versus CSE",
  CRP     = "CRP versus CSP",
  VRE     = "VRE versus VSE",
  MRSA    = "MRSA versus MSSA"
)

# --------------------------------------------------
# Clean data
# --------------------------------------------------
df_count_all <- df_count_all %>%
  filter(pathogen %in% keep_pathogen) %>%
  mutate(
    pathogen_label = recode(pathogen, !!!pathogen_map),
    resistance = case_when(
      grepl("sus", group, ignore.case = TRUE) ~ "S",
      grepl("res", group, ignore.case = TRUE) ~ "R",
      TRUE ~ NA_character_
    ),
    level_clean = case_when(
      level %in% c("overall", "Overall") ~ "Overall",
      
      level %in% c("VAP") ~ "VAP",
      level %in% c("Hospital-acquired BSI", "Hospital−acquired BSI") ~ "Hospital-acquired BSI",
      level %in% c("Healthcare-associated BSI", "Healthcare−associated BSI") ~ "Healthcare-associated BSI",
      
      level %in% c("<1 year", "< 1 year", "Infants aged <1 year") ~ "<1 year",
      level %in% c("1–4 years", "1-4 years", "1 to 4 years") ~ "1–4 years",
      level %in% c("5–14 years", "5-14 years", "5 to 14 years") ~ "12–14 years",
      level %in% c("15–49 years", "15-49 years", "15 to 49 years") ~ "15–49 years",
      level %in% c("50–69 years", "50-69 years", "50 to 69 years") ~ "50–69 years",
      level %in% c("≥70 years", ">=70 years", "70+ years", "≥ 70 years") ~ "≥70 years",
      
      level %in% c("High income") ~ "High income",
      level %in% c("Upper middle income") ~ "Upper middle income",
      level %in% c("Lower middle income") ~ "Lower middle income",
      
      TRUE ~ as.character(level)
    ),
    section = case_when(
      subgroup == "overall"         ~ "Overall",
      subgroup == "infection_types" ~ "Infection syndromes",
      subgroup == "age_group_new"   ~ "Age group",
      subgroup == "country_income"  ~ "World Bank income status",
      TRUE ~ subgroup
    ),
    section_order = case_when(
      subgroup == "overall"         ~ 1,
      subgroup == "infection_types" ~ 2,
      subgroup == "age_group_new"   ~ 3,
      subgroup == "country_income"  ~ 4,
      TRUE ~ 99
    ),
    within_order = case_when(
      subgroup == "overall" ~ 1,
      
      subgroup == "infection_types" & level_clean == "VAP" ~ 1,
      subgroup == "infection_types" & level_clean == "Hospital-acquired BSI" ~ 2,
      subgroup == "infection_types" & level_clean == "Healthcare-associated BSI" ~ 3,
      
      subgroup == "age_group_new" & level_clean == "<1 year" ~ 1,
      subgroup == "age_group_new" & level_clean == "1–4 years" ~ 2,
      subgroup == "age_group_new" & level_clean == "12–14 years" ~ 3,
      subgroup == "age_group_new" & level_clean == "15–49 years" ~ 4,
      subgroup == "age_group_new" & level_clean == "50–69 years" ~ 5,
      subgroup == "age_group_new" & level_clean == "≥70 years" ~ 6,
      
      subgroup == "country_income" & level_clean == "High income" ~ 1,
      subgroup == "country_income" & level_clean == "Upper middle income" ~ 2,
      subgroup == "country_income" & level_clean == "Lower middle income" ~ 3,
      
      TRUE ~ 99
    )
  )

# --------------------------------------------------
# Build display table: ONE display_order per subgroup line
# --------------------------------------------------
build_eq_table_rs <- function(df, pathogen_name) {
  
  dat <- df %>%
    filter(pathogen_label == pathogen_name) %>%
    filter(!is.na(resistance)) %>%
    arrange(section_order, within_order, resistance)
  
  out_list <- list()
  line_id <- 1
  
  # overall
  dat_overall <- dat %>% filter(section == "Overall")
  if (nrow(dat_overall) > 0) {
    out_list[[length(out_list) + 1]] <- dat_overall %>%
      transmute(
        pathogen_label,
        subgroup_label = "Overall",
        resistance,
        adjusted_mean,
        lower_ci,
        upper_ci,
        is_header = FALSE,
        display_order = line_id
      )
    line_id <- line_id + 1
  }
  
  # infection header
  dat_inf <- dat %>% filter(section == "Infection syndromes")
  if (nrow(dat_inf) > 0) {
    out_list[[length(out_list) + 1]] <- data.frame(
      pathogen_label = pathogen_name,
      subgroup_label = "Infection syndromes",
      resistance = NA_character_,
      adjusted_mean = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      is_header = TRUE,
      display_order = line_id,
      stringsAsFactors = FALSE
    )
    line_id <- line_id + 1
    
    for (lv in unique(dat_inf$level_clean)) {
      tmp <- dat_inf %>% filter(level_clean == lv)
      out_list[[length(out_list) + 1]] <- tmp %>%
        transmute(
          pathogen_label,
          subgroup_label = paste0("   ", lv),
          resistance,
          adjusted_mean,
          lower_ci,
          upper_ci,
          is_header = FALSE,
          display_order = line_id
        )
      line_id <- line_id + 1
    }
  }
  
  # age header
  dat_age <- dat %>% filter(section == "Age group")
  if (nrow(dat_age) > 0) {
    out_list[[length(out_list) + 1]] <- data.frame(
      pathogen_label = pathogen_name,
      subgroup_label = "Age group",
      resistance = NA_character_,
      adjusted_mean = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      is_header = TRUE,
      display_order = line_id,
      stringsAsFactors = FALSE
    )
    line_id <- line_id + 1
    
    for (lv in unique(dat_age$level_clean)) {
      tmp <- dat_age %>% filter(level_clean == lv)
      out_list[[length(out_list) + 1]] <- tmp %>%
        transmute(
          pathogen_label,
          subgroup_label = paste0("   ", lv),
          resistance,
          adjusted_mean,
          lower_ci,
          upper_ci,
          is_header = FALSE,
          display_order = line_id
        )
      line_id <- line_id + 1
    }
  }
  
  # income header
  dat_inc <- dat %>% filter(section == "World Bank income status")
  if (nrow(dat_inc) > 0) {
    out_list[[length(out_list) + 1]] <- data.frame(
      pathogen_label = pathogen_name,
      subgroup_label = "World Bank income status",
      resistance = NA_character_,
      adjusted_mean = NA_real_,
      lower_ci = NA_real_,
      upper_ci = NA_real_,
      is_header = TRUE,
      display_order = line_id,
      stringsAsFactors = FALSE
    )
    line_id <- line_id + 1
    
    for (lv in unique(dat_inc$level_clean)) {
      tmp <- dat_inc %>% filter(level_clean == lv)
      out_list[[length(out_list) + 1]] <- tmp %>%
        transmute(
          pathogen_label,
          subgroup_label = paste0("   ", lv),
          resistance,
          adjusted_mean,
          lower_ci,
          upper_ci,
          is_header = FALSE,
          display_order = line_id
        )
      line_id <- line_id + 1
    }
  }
  
  bind_rows(out_list)
}

# --------------------------------------------------
# Build plot base
# --------------------------------------------------
pathogen_order <- c("CRA", "3GCRE", "CRE", "CRP", "VRE", "MRSA")

plot_base <- bind_rows(
  lapply(pathogen_order, function(p) build_eq_table_rs(df_count_all, p))
)

# --------------------------------------------------
# Shared x limits
# --------------------------------------------------
get_xlim <- function(dat) {
  xmin <- min(dat$lower_ci, na.rm = TRUE)
  xmax <- max(dat$upper_ci, na.rm = TRUE)
  pad  <- 0.08 * (xmax - xmin)
  c(xmin - pad, xmax + pad)
}

xlim_all <- get_xlim(plot_base)

# --------------------------------------------------
# Plot function
# --------------------------------------------------
make_eq_plot_rs <- function(dat, pathogen_name, xlim_use,
                            show_group_labels = TRUE,
                            show_x_title = FALSE) {
  
  dfp <- dat %>%
    filter(pathogen_label == pathogen_name) %>%
    arrange(display_order, resistance)
  
  y_mult <- 1.28
  label_size <- 4.0
  header_size <- 4.1
  point_size <- 1.8
  title_size <- 12
  axis_text_size <- 11
  
  y_map <- dfp %>%
    distinct(display_order, subgroup_label, is_header) %>%
    arrange(display_order) %>%
    mutate(y = -row_number() * y_mult)
  
  dfp <- dfp %>%
    left_join(y_map, by = c("display_order", "subgroup_label", "is_header")) %>%
    mutate(
      y_plot = case_when(
        is_header ~ y,
        resistance == "R" ~ y + 0.08,
        resistance == "S" ~ y - 0.08,
        TRUE ~ y
      ),
      is_bold = subgroup_label %in% c(
        "Overall",
        "Infection syndromes",
        "Age group",
        "World Bank income status"
      )
    )
  
  section_lines <- y_map %>%
    filter(subgroup_label %in% c(
      "Infection syndromes",
      "Age group",
      "World Bank income status"
    ))
  
  subgroup_x <- xlim_use[1] - 0.14 * diff(xlim_use)
  left_margin <- if (show_group_labels) 190 else 30
  header_y <- max(y_map$y) + 1.5
  
  p <- ggplot() +
    geom_hline(
      data = section_lines,
      aes(yintercept = y),
      colour = "grey82",
      linewidth = 0.25
    ) +
    # horizontal CI line
    geom_segment(
      data = dfp %>% filter(!is_header),
      aes(
        x = lower_ci, xend = upper_ci,
        y = y_plot,  yend = y_plot,
        color = resistance
      ),
      linewidth = 0.55
    ) +
    # left cap
    geom_segment(
      data = dfp %>% filter(!is_header),
      aes(
        x = lower_ci, xend = lower_ci,
        y = y_plot - 0.18, yend = y_plot + 0.18,
        color = resistance
      ),
      linewidth = 0.55
    ) +
    # right cap
    geom_segment(
      data = dfp %>% filter(!is_header),
      aes(
        x = upper_ci, xend = upper_ci,
        y = y_plot - 0.18, yend = y_plot + 0.18,
        color = resistance
      ),
      linewidth = 0.55
    ) +
    geom_point(
      data = dfp %>% filter(!is_header),
      aes(x = adjusted_mean, y = y_plot, color = resistance),
      size = point_size
    ) +
    coord_cartesian(
      xlim = xlim_use,
      ylim = c(min(y_map$y) - 0.35, header_y - 0.8),
      clip = "off"
    ) +
    scale_color_manual(
      values = c("R" = "#1f77b4", "S" = "#E3B505"),
      breaks = c("R", "S"),
      labels = c("Resistant (R)", "Susceptible (S)")
    ) +
    labs(
      title = title_map[[pathogen_name]],
      x = if (show_x_title) "Adjusted mean EQ-5D-3L utility values" else NULL,
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 12, base_family = "Times") +
    theme(
      plot.title = element_text(
        face = "bold",
        size = title_size,
        hjust = 0.5,
        margin = margin(b = 18)
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(
        size = axis_text_size,
        color = "black",
        margin = margin(t = 25)
      ),
      axis.text.x = element_text(
        size = axis_text_size,
        color = "black"
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey85"),
      plot.margin = margin(8, 8, 10, left_margin),
      legend.text = element_text(size = 12, family = "Times"),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "bottom"
    )
  
  # arrow + labels under x-axis
  arrow_y <- min(y_map$y) - 2.5
  
  p <- p +
    annotate(
      "segment",
      x = xlim_use[1]-0.1,
      xend = xlim_use[2]+0.12,
      y = arrow_y,
      yend = arrow_y,
      arrow = arrow(ends = "both", type = "closed", length = unit(0.15, "cm")),
      linewidth = 0.4,
      color = "grey40"
    ) +
    annotate(
      "text",
      x = xlim_use[1]-0.1,
      y = arrow_y - 0.55,
      label = "Worse health",
      hjust = 0,
      size = 3.5,
      family = "Times",
      color = "grey40"
    ) +
    annotate(
      "text",
      x = xlim_use[2]+0.12,
      y = arrow_y - 0.55,
      label = "Better health",
      hjust = 1,
      size = 3.5,
      family = "Times",
      color = "grey40"
    )
  
  
  if (show_group_labels) {
    label_df <- y_map %>%
      mutate(
        fontface_lab = ifelse(
          subgroup_label %in% c("Overall", "Infection syndromes", "Age group", "World Bank income status"),
          "bold", "plain"
        )
      )
    
    p <- p +
      geom_text(
        data = label_df,
        aes(x = subgroup_x, y = y, label = subgroup_label, fontface = fontface_lab),
        hjust = 1,
        size = label_size,
        family = "Times",
        inherit.aes = FALSE
      ) +
      annotate(
        "text",
        x = subgroup_x,
        y = header_y,
        label = "Subgroups",
        hjust = 1,
        fontface = "bold",
        family = "Times",
        size = header_size
      )
  }
  
  p
}
# --------------------------------------------------
# Build 6 panels
# --------------------------------------------------
fig_eq_rs1 <- wrap_plots(
  list(
    make_eq_plot_rs(plot_base, "CRA",   xlim_all, TRUE,  FALSE),
    make_eq_plot_rs(plot_base, "3GCRE", xlim_all, FALSE, TRUE),
    make_eq_plot_rs(plot_base, "CRE",   xlim_all, FALSE, FALSE)
  ),
  ncol = 3,
  guides = "collect"
) & theme(
  legend.position = "bottom",
  legend.direction = "horizontal"
)

fig_eq_rs2 <- wrap_plots(
  list(
    make_eq_plot_rs(plot_base, "CRP",   xlim_all, TRUE,  FALSE),
    make_eq_plot_rs(plot_base, "VRE",   xlim_all, FALSE, TRUE),
    make_eq_plot_rs(plot_base, "MRSA",  xlim_all, FALSE, FALSE)
  ),
  ncol = 3,
  guides = "collect"
) & theme(
  legend.position = "bottom",
  legend.direction = "horizontal"
)

# --------------------------------------------------
# Save
# --------------------------------------------------
CairoPDF(
  file = "output/eq5d3l_RS1.pdf",
  width = 13,
  height = 6
)
print(fig_eq_rs1)
dev.off()

CairoPDF(
  file = "output/eq5d3l_RS2.pdf",
  width = 13,
  height = 6
)
print(fig_eq_rs2)
dev.off()

