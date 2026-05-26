# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 survival,
                 survminer,
                 ggplot2,
                 gridExtra,
                 ggpubr,
                 extrafont,
                 patchwork,
                 Cairo,
                 stringr)
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
         "pse_car", "entc_van", "sa_meth")

# Raw labels for group assignment
raw_labels <- list(
  c("CSA", "CRA"),
  c("3GCSE", "3GCRE"),
  c("CSE", "CRE"),
  c("CSP", "CRP"),
  c("VSE", "VRE"),
  c("MSSA", "MRSA")
)

# Create padded labels for visual alignment
max_len <- max(nchar(unlist(raw_labels)))
labels <- lapply(raw_labels, function(x) str_pad(x, max_len, side = "left"))

# Initialize list to store plots
combined_plot <- list()

# Loop over each dataset
for (i in 1:6) {
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

  # Fit survival curves
  fit <- survfit(Surv(time, event) ~ country_income + ris_group, data = df)
  
  # Calculate log-rank p-values by income group
  df_split <- df %>% group_split(country_income)
  p_vals <- lapply(df_split, function(sub_df) {
    survdiff(Surv(time, event) ~ ris_group, data = sub_df)
  }) %>%
    lapply(function(p_model) {
      df_non_empty <- sum(p_model$n > 0) - 1
      1 - pchisq(p_model$chisq, df = df_non_empty)
    })
  
  p_labels <- sapply(p_vals, function(pv) {
    if (is.na(pv)) return("")
    
    if (pv < 0.0001) {
      "< 0.0001"
    } else {
      paste0("= ", sprintf("%.4f", pv))
    }
  })
  
  pvals_df <- data.frame(
    country_income = c("High income", "Upper middle income", "Lower middle income"),
    x = c(4, 4, 4),
    y = c(0.2, 0.2, 0.2),
    label = paste0("p ", p_labels)
  )
  
  # Get survival object with risk table
  ggsurv <- ggsurvplot(fit, 
                       data = df,
                       risk.table = TRUE, 
                       break.time.by = 7,
                       xlim = c(0, 28),
                       tables.theme = theme_minimal())
  
  # Extract risk table data and fix group labeling
  tbl_data <- ggsurv$table$data %>%
    mutate(
      ris_group = sub(".*, ris_group=", "", strata),
      ris_group = sub("\\)", "", ris_group),
      ris_label = factor(dplyr::recode(ris_group, !!!setNames(labels[[i]], raw_labels[[i]])),
                         levels = labels[[i]])
    ) %>%
    filter(!is.na(time), !is.na(n.risk), !is.na(ris_label), time <= 28)
  
  # Plot risk table
  tbl_fct <- ggplot(tbl_data, aes(time, ris_label, color = ris_label)) + 
    geom_text(aes(label = n.risk), size = 6, family = "Times New Roman") +
    facet_wrap(~country_income, labeller = label_both) +  
    theme_minimal(base_family = "Times New Roman") + 
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5), 
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
      strip.background = element_rect(fill = "white", linewidth = 1),
      strip.text = element_blank(),
      panel.spacing = unit(0.6, "cm"),
      axis.text.y = element_text(color = "black", size = 16, 
                                 family = "Times New Roman", 
                                 margin = margin(r = ifelse(i == 4, 20, 16))),
      axis.text.x = element_blank(),
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
    ggtitle(" ") +
    xlab("") + 
    ylab("")
  
  # Plot survival curves
  # Modify country_income order with prefix
  income_levels <- c("1. High income", "2. Upper middle income", "3. Lower middle income")
  names(income_levels) <- c("High income", "Upper middle income", "Lower middle income")
  df$country_income <- factor(income_levels[as.character(df$country_income)],
                              levels = income_levels)
  
  # p-value text coordinate
  pvals_df$country_income <- factor(income_levels[as.character(pvals_df$country_income)],
                                    levels = income_levels)
  
  # Plot
  plt_fct <- ggsurvplot_facet(
    fit,
    data = df,
    facet.by = "country_income",
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
        axis.text.y = element_text(
          color = "black",
          size = 16,
          family = "Times New Roman",
          hjust = 1,
          margin = margin(r = 14)),
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
    panel.labs = list(country_income = income_levels), 
    panel.labs.font = list(size = 16, family = "Times New Roman"),
    short.panel.labs = TRUE,
    legend.labs = labels[[i]],
    legend = "none"
  ) +
    geom_text(data = pvals_df, aes(x = x, y = y, label = label),
              size = 6, family = "Times New Roman") +
    facet_wrap(~country_income, labeller = labeller(country_income = function(x) sub("^[0-9]+\\. ", "", x)))
  
  # Combine survival plot and table
  combined_plot[[i]] <- (plt_fct + theme(plot.margin = margin(t = 2, r = 1, b = 2, l = 1))) /
    (tbl_fct + theme(plot.margin = margin(t = 0, r = 1, b = 2, l = 1))) +
    plot_layout(heights = c(7.5, 2.5)) & 
    theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))
  
  print(combined_plot[[i]])
}

# Save all plots
pdf_files <- paste0("output/sub_figure/att28_death_KM_0", 1:6, "_income.pdf")

for (i in seq_along(combined_plot)) {
  CairoPDF(file = pdf_files[i], width = 15, height = 6)
  print(combined_plot[[i]])
  dev.off()
}

