# Clear environment
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr, magrittr, survival, survminer, ggplot2,
    gridExtra, ggpubr, extrafont, patchwork, Cairo,
    stringr, grid, gtable, ggplotify
  )
})

# Working directory & data
wd <- "./"
setwd(wd)

df_all <- readRDS("data/att_first28.RData")

# Fonts
loadfonts()

# Variables & labels
var <- c("aci_car", "ent_thir", "ent_car",
         "pse_car", "entc_van", "sa_meth")

raw_labels <- list(
  c("CSA", "CRA"),
  c("3GCSE", "3GCRE"),
  c("CSE", "CRE"),
  c("CSP", "CRP"),
  c("VSE", "VRE"),
  c("MSSA", "MRSA")
)

# Pad labels with leading spaces for alignment
max_len <- max(nchar(unlist(raw_labels)))
labels <- lapply(raw_labels, function(x) stringr::str_pad(x, max_len, side = "left"))

# Panel order
full_order <- c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")

# A “truly empty” panel (white background, no elements),
# used as a placeholder for VRE
empty_panel <- ggplot() +
  theme_void() +
  ggtitle(" ") +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.background = element_rect(fill = "white", color = NA)
  )


# Main loop
combined_plot <- list()

for (i in 1:6) {
  df <- df_all[[i]]
  ris_name <- var[i]
  
  # --------
  # For VRE: do not forcibly add VAP; set reference as H-acq BSI
  # --------
  if (ris_name == "entc_van") {
    df$infection_types <- droplevels(factor(df$infection_types))
    df$infection_types <- relevel(df$infection_types, ref = "Hospital-acquired BSI")
  }
  
  # Resistant vs susceptible grouping and labels
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
  
  # Shift 28-day events to 27.999 to avoid boundary issue
  df$time[df$time == 28 & df$event == 1] <- 27.999
  
  # --------
  # Panel levels and labels
  # For non-VRE: typically all three panels
  # For VRE: only use panels actually present
  # --------
  present_levels <- levels(droplevels(factor(df$infection_types)))
  
  if (ris_name == "entc_van") {
    use_levels <- intersect(full_order, present_levels)
  } else {
    use_levels <- intersect(full_order, present_levels)
  }
  
  # Add numeric prefix to labels: e.g., "1. VAP"
  infection_levels <- setNames(paste0(seq_along(use_levels), ". ", use_levels), use_levels)
  df$infection_types <- factor(infection_levels[as.character(df$infection_types)],
                               levels = unname(infection_levels))
  
  # --------
  # Fit KM curves
  # --------
  fit <- survfit(Surv(time, event) ~ infection_types + ris_group, data = df)
  
  # Log-rank p-values per panel
  df_split <- df %>% dplyr::group_split(infection_types)
  p_vals <- lapply(df_split, function(sub_df) {
    if (nrow(sub_df) == 0 || all(is.na(sub_df$ris_group))) return(NA_real_)
    pm <- survdiff(Surv(time, event) ~ ris_group, data = sub_df)
    df_non_empty <- sum(pm$n > 0) - 1
    if (df_non_empty < 1) return(NA_real_)
    1 - pchisq(pm$chisq, df = df_non_empty)
  })
  
  p_labels <- sapply(p_vals, function(pv) {
    if (is.na(pv)) return("")
    
    if (pv < 0.0001) {
      "< 0.0001"
    } else {
      paste0("= ", sprintf("%.4f", pv))
    }
  })
  
  split_lvls <- sapply(df_split, function(d) sub("^\\d+\\.\\s*", "", as.character(unique(d$infection_types))))
  names(p_labels) <- split_lvls
  
  show_levels <- names(infection_levels)[names(infection_levels) %in% names(p_labels) & p_labels[names(infection_levels)] != ""]
  pvals_df <- data.frame(
    infection_types = unname(infection_levels[show_levels]),
    x = 4, y = 0.2,
    label = paste0("p ", p_labels[show_levels]),
    check.names = FALSE
  )
  
  # =========================
  # Risk table
  # =========================
  ggsurv <- ggsurvplot(fit,
                       data = df,
                       risk.table = TRUE,
                       break.time.by = 7,
                       xlim = c(0, 28),
                       tables.theme = theme_minimal())
  
  tbl_data <- ggsurv$table$data %>%
    mutate(
      ris_group = sub(".*, ris_group=", "", strata),
      ris_group = sub("\\)", "", ris_group),
      ris_label = factor(dplyr::recode(ris_group, !!!setNames(labels[[i]], raw_labels[[i]])),
                         levels = labels[[i]]),
      infection_types = factor(as.character(.data$infection_types),
                               levels = unname(infection_levels))
    ) %>%
    filter(!is.na(time), !is.na(n.risk), !is.na(ris_label), time <= 28)
  
  tbl_plot <- ggplot(tbl_data, aes(time, ris_label, color = ris_label)) +
    geom_text(aes(label = n.risk), size = 6, family = "Times New Roman") +
    facet_wrap(~infection_types,
               labeller = labeller(infection_types = function(x) sub("^[0-9]+\\. ", "", x)),
               drop = FALSE, nrow = 1) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
      strip.background = element_rect(fill = "white", linewidth = 1),
      strip.text = element_blank(),
      panel.spacing = unit(0.6, "cm"),
      axis.text.y = element_text(color = "black", size = 16,
                                 family = "Times New Roman",
                                 margin = margin(r = ifelse(i == 2, 12, 15))),
      axis.text.x = element_text(color = "black", size = 14,
                                 family = "Times New Roman"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      legend.position = "none"
    ) +
    scale_y_discrete(limits = labels[[i]]) +
    scale_x_continuous(
      limits = c(0, 28),
      breaks = c(0, 7, 14, 21, 28),
      labels = c("0", "7", "14", "21", "28")
    ) +
    scale_colour_manual(values = setNames(c("#EFC000FF", "#0073c2FF"), labels[[i]])) +
    ggtitle(" ") + xlab("") + ylab("")
  
  # =========================
  # Survival curves
  # =========================
  surv_plot_base <- ggsurvplot_facet(
    fit,
    data = df,
    facet.by = "infection_types",
    pval = FALSE,
    palette = c("#EFC000FF", "#0073c2FF"),
    conf.int = TRUE,
    linetype = 1,
    ggtheme = theme_minimal(base_family = "Times New Roman") +
      theme(
        panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
        panel.spacing = unit(0.6, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.1, "cm"),
        axis.text.x = element_text(color = "black", size = 16, family = "Times New Roman"),
        axis.text.y = element_text(color = "black", margin = margin(r = 12),
                                   size = 16, family = "Times New Roman"),
        axis.title.x = element_text(color = "black", margin = margin(t = 15), size = 16,
                                    family = "Times New Roman"),
        axis.title.y = element_text(color = "black", margin = margin(r = 12), size = 16,
                                    family = "Times New Roman"),
        strip.text.x = element_text(size = 16, margin = margin(t = 6, b = 6)),
        strip.background.x = element_rect(color = "black", fill = "white", linewidth = NA)
      ),
    title = " ",
    xlab = "Follow-up time since infection onset (days)",
    ylab = "Survival probability",
    font.tickslab = c(16, "plain"),
    break.time.by = 7,
    xlim = c(0, 28),
    panel.labs = list(infection_types = infection_levels),
    panel.labs.font = list(size = 16, family = "Times New Roman"),
    short.panel.labs = TRUE,
    legend.labs = labels[[i]],
    legend = "none"
  ) +
    # Add panel-wise p-values (only where data exist)
    geom_text(data = pvals_df, aes(x = x, y = y, label = label),
              size = 6, family = "Times New Roman") +
    facet_wrap(
      ~infection_types,
      labeller = labeller(infection_types = function(x) sub("^[0-9]+\\. ", "", x)),
      drop = FALSE,
      nrow = 1
    )
  
  # Convert to ggplot objects (for patchwork composition)
  top_plot    <- ggplotify::as.ggplot(surv_plot_base)
  bottom_plot <- ggplotify::as.ggplot(tbl_plot)
  
  # =========================
  # Combine plots
  # Non-VRE: normal layout
  # VRE: add blank panel on the right
  # =========================
  if (ris_name == "entc_van") {
    # Top row: KM + blank panel
    top_row <- (top_plot | empty_panel) + plot_layout(widths = c(2, 1))
    # Bottom row: risk table + blank panel
    bottom_row <- (bottom_plot | empty_panel) + plot_layout(widths = c(2, 1))
    
    combined_plot[[i]] <- (top_row + theme(plot.margin = margin(t = 2, r = 1, b = 2, l = 1))) /
      (bottom_row + theme(plot.margin = margin(t = 0, r = 1, b = 2, l = 1))) +
      plot_layout(heights = c(7.5, 2.5)) &
      theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))
    
  } else {
    combined_plot[[i]] <- (top_plot + theme(plot.margin = margin(t = 2, r = 1, b = 2, l = 1))) /
      (bottom_plot + theme(plot.margin = margin(t = 0, r = 1, b = 2, l = 1))) +
      plot_layout(heights = c(7.5, 2.5)) &
      theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))
  }
  
  print(combined_plot[[i]])
}


# Save as PDF
pdf_files <- paste0("output/sub_figure/att28_death_KM_0", 1:6, "_infection.pdf")
for (i in seq_along(combined_plot)) {
  CairoPDF(file = pdf_files[i], width = 15, height = 6)
  print(combined_plot[[i]])
  dev.off()
}

