# Clear environment
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr, tidyr, magrittr, survival, survminer, ggplot2,
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

# Panel order for age groups
full_order <- c(
  "<1 year",
  "1–4 years",
  "5–14 years",
  "15–49 years",
  "50–69 years",
  "≥70 years"
)

# color map
get_col_map <- function(i) {
  setNames(c("#EFC000FF", "#0073c2FF"), raw_labels[[i]])
}

# p-value formatter
fmt_p <- function(pv) {
  if (is.na(pv)) return("")
  
  if (pv < 0.0001) {
    return("p < 0.0001")
  }
  
  paste0("p = ", sprintf("%.4f", pv))
}

combined_plot <- vector("list", 6)

for (i in 1:6) {
  df <- df_all[[i]]
  ris_name <- var[i]
  
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
  
  # Keep age group order
  df$age_group_new <- factor(as.character(df$age_group_new), levels = full_order)
  
  # Keep only age groups where both groups > 10
  keep_age <- df %>%
    filter(!is.na(age_group_new), !is.na(ris_group)) %>%
    count(age_group_new, ris_group) %>%
    pivot_wider(names_from = ris_group, values_from = n, values_fill = 0)
  
  keep_levels <- keep_age %>%
    filter(
      .data[[raw_labels[[i]][1]]] > 10,
      .data[[raw_labels[[i]][2]]] > 10
    ) %>%
    pull(age_group_new) %>%
    as.character()
  
  df <- df %>%
    filter(age_group_new %in% keep_levels) %>%
    filter(!is.na(age_group_new), !is.na(ris_group), !is.na(time), !is.na(event))
  
  if (nrow(df) == 0) {
    combined_plot[[i]] <- ggplot() + theme_void()
    print(combined_plot[[i]])
    next
  }
  
  # Shift 28-day events to 27.999 to avoid boundary issue
  df$time[df$time == 28 & df$event == 1] <- 27.999
  
  # =========================================================
  # FIGURE 1-3: each curve directly above its own table
  # =========================================================
  if (i <= 3) {
    
    age_present <- intersect(full_order, unique(as.character(df$age_group_new)))
    panel_list <- list()
    
    for (j in seq_along(age_present)) {
      ag <- age_present[j]
      show_y <- ((j - 1) %% 3 == 0)
      
      sub_df <- df %>%
        filter(age_group_new == ag) %>%
        droplevels()
      
      if (n_distinct(sub_df$ris_group) < 2) next
      
      fit_sub <- survfit(Surv(time, event) ~ ris_group, data = sub_df)
      
      pv <- tryCatch({
        sdiff <- survdiff(Surv(time, event) ~ ris_group, data = sub_df)
        df_non_empty <- sum(sdiff$n > 0) - 1
        if (df_non_empty < 1) NA_real_ else 1 - pchisq(sdiff$chisq, df = df_non_empty)
      }, error = function(e) NA_real_)
      
      p_lab <- fmt_p(pv)
      
      g_sub <- ggsurvplot(
        fit_sub,
        data = sub_df,
        risk.table = TRUE,
        conf.int = TRUE,
        pval = FALSE,
        break.time.by = 7,
        xlim = c(0, 28),
        legend = "none",
        palette = unname(get_col_map(i)),
        linetype = 1,
        censor = FALSE,
        ggtheme = theme_minimal(base_family = "Times New Roman"),
        tables.theme = theme_minimal(base_family = "Times New Roman")
      )
      
      # top curve
      top_plot <- g_sub$plot +
        labs(
          title = ag,
          x = "Follow-up time since infection onset (days)",
          y = if (show_y) "Survival probability" else NULL
        ) +
        annotate("text", x = 4, y = 0.2, label = p_lab,
                 family = "Times New Roman", size = 4.5, color = "black") +
        theme(
          panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(color = "black", size = 14, family = "Times New Roman"),
          axis.text.y = if (show_y) {
            element_text(color = "black", size = 14, family = "Times New Roman", margin = margin(r = 4))
          } else {
            element_blank()
          },
          axis.title.x = element_text(color = "black", size = 14, family = "Times New Roman",
                                      margin = margin(t = 10)),
          axis.title.y = if (show_y) {
            element_text(color = "black", size = 14, family = "Times New Roman",
                         margin = margin(r = 4))
          } else {
            element_blank()
          },
          axis.ticks.y = if (show_y) {
            element_line(color = "black", linewidth = 0.5)
          } else {
            element_blank()
          },
          plot.title = element_text(
            size = 15, color = "black", family = "Times New Roman",
            hjust = 0.5, margin = margin(t = 6, b = 6)
          ),
          plot.title.position = "plot",
          plot.margin = margin(t = 4, r = 4, b = 2, l = 4)
        )
      
      # bottom risk table
      tbl_data <- g_sub$table$data %>%
        mutate(
          ris_group = sub(".*=", "", strata),
          ris_group = factor(ris_group, levels = raw_labels[[i]]),
          ris_label = factor(
            dplyr::recode(as.character(ris_group), !!!setNames(labels[[i]], raw_labels[[i]])),
            levels = labels[[i]]
          )
        ) %>%
        filter(time <= 28)
      
      bottom_plot <- ggplot(tbl_data, aes(time, ris_label, color = ris_group)) +
        geom_text(aes(label = n.risk), size = 5, family = "Times New Roman") +
        theme_minimal(base_family = "Times New Roman") +
        theme(
          panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
          axis.text.y = if (show_y) {
            element_text(color = "black", size = 14, family = "Times New Roman",
                         margin = margin(r = 12))
          } else {
            element_blank()
          },
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = if (show_y) {
            element_line(color = "black", linewidth = 0.5)
          } else {
            element_blank()
          },
          axis.ticks.length = unit(0.1, "cm"),
          legend.position = "none",
          plot.margin = margin(t = 0, r = 4, b = 2, l = 4)
        ) +
        scale_y_discrete(limits = labels[[i]]) +
        scale_x_continuous(
          limits = c(0, 28),
          breaks = c(0, 7, 14, 21, 28),
          labels = rep("", 5)
        ) +
        scale_colour_manual(values = get_col_map(i)) +
        coord_cartesian(clip = "off")
      
      panel_list[[ag]] <- top_plot / bottom_plot +
        plot_layout(heights = c(7.5, 2.5))
    }
    
    panel_list <- panel_list[age_present]
    
    n_missing <- (ceiling(length(panel_list) / 3) * 3) - length(panel_list)
    if (n_missing > 0) {
      for (k in seq_len(n_missing)) {
        panel_list[[paste0("blank_", k)]] <- patchwork::plot_spacer()
      }
    }
    
    combined_plot[[i]] <- wrap_plots(panel_list, ncol = 3) &
      theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))
    
    print(combined_plot[[i]])
  }
  
  # =========================================================
  # FIGURE 4-6:
  # single panels + manual blank panel for Figure 5
  # =========================================================
  if (i > 3) {
    
    age_present <- intersect(full_order, unique(as.character(df$age_group_new)))
    top_panels <- list()
    bottom_panels <- list()
    panel_idx <- 0
    
    for (ag in age_present) {
      
      sub_df <- df %>%
        filter(age_group_new == ag) %>%
        droplevels()
      
      if (nrow(sub_df) == 0 || dplyr::n_distinct(sub_df$ris_group) < 2) next
      
      panel_idx <- panel_idx + 1
      show_y <- (panel_idx == 1)
      
      fit_sub <- survfit(Surv(time, event) ~ ris_group, data = sub_df)
      
      pv <- tryCatch({
        sdiff <- survdiff(Surv(time, event) ~ ris_group, data = sub_df)
        df_non_empty <- sum(sdiff$n > 0) - 1
        if (df_non_empty < 1) NA_real_ else 1 - pchisq(sdiff$chisq, df = df_non_empty)
      }, error = function(e) NA_real_)
      
      p_lab <- fmt_p(pv)
      
      g_sub <- ggsurvplot(
        fit_sub,
        data = sub_df,
        risk.table = TRUE,
        conf.int = TRUE,
        pval = FALSE,
        break.time.by = 7,
        xlim = c(0, 28),
        legend = "none",
        palette = unname(get_col_map(i)),
        linetype = 1,
        censor = FALSE,
        ggtheme = theme_minimal(base_family = "Times New Roman"),
        tables.theme = theme_minimal(base_family = "Times New Roman")
      )
      
      # ---------- top KM ----------
      top_plot_sub <- g_sub$plot +
        labs(
          title = ag,
          x = "Follow-up time since infection onset (days)",
          y = if (show_y) "Survival probability" else NULL
        ) +
        annotate(
          "text",
          x = 4, y = 0.2,
          label = p_lab,
          family = "Times New Roman",
          size = 4.5,
          color = "black"
        ) +
        theme(
          panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length = unit(0.1, "cm"),
          axis.text.x = element_text(color = "black", size = 16, family = "Times New Roman"),
          axis.text.y = if (show_y) {
            element_text(color = "black", size = 16, family = "Times New Roman", margin = margin(r = 12))
          } else {
            element_blank()
          },
          axis.title.x = element_text(
            color = "black", size = 16, family = "Times New Roman",
            margin = margin(t = 15)
          ),
          axis.title.y = if (show_y) {
            element_text(color = "black", size = 16, family = "Times New Roman",
                         margin = margin(r = 8))
          } else {
            element_blank()
          },
          axis.ticks.y = if (show_y) {
            element_line(color = "black", linewidth = 0.5)
          } else {
            element_blank()
          },
          plot.title = element_text(
            size = 16, color = "black", family = "Times New Roman",
            hjust = 0.5, margin = margin(t = 6, b = 6)
          ),
          plot.title.position = "plot",
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
        )
      
      # ---------- bottom risk table ----------
      tbl_data_sub <- g_sub$table$data %>%
        mutate(
          ris_group = sub(".*=", "", strata),
          ris_group = factor(ris_group, levels = raw_labels[[i]]),
          ris_label = factor(
            dplyr::recode(as.character(ris_group), !!!setNames(labels[[i]], raw_labels[[i]])),
            levels = labels[[i]]
          )
        ) %>%
        filter(time <= 28)
      
      bottom_plot_sub <- ggplot(tbl_data_sub, aes(time, ris_label, color = ris_group)) +
        geom_text(aes(label = n.risk), size = 6, family = "Times New Roman") +
        theme_minimal(base_family = "Times New Roman") +
        theme(
          panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
          axis.text.y = if (show_y) {
            element_text(
              color = "black",
              size = 16,
              family = "Times New Roman",
              margin = margin(r = 14)
            )
          } else {
            element_blank()
          },
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = if (show_y) {
            element_line(color = "black", linewidth = 0.5)
          } else {
            element_blank()
          },
          axis.ticks.length = unit(0.1, "cm"),
          legend.position = "none",
          plot.margin = margin(t = 0, r = 2, b = 2, l = 2)
        ) +
        scale_y_discrete(limits = labels[[i]]) +
        scale_x_continuous(
          limits = c(0, 28),
          breaks = c(0, 7, 14, 21, 28),
          labels = rep("", 5)
        ) +
        scale_colour_manual(values = get_col_map(i)) +
        coord_cartesian(clip = "off")
      
      top_panels[[ag]] <- top_plot_sub
      bottom_panels[[ag]] <- bottom_plot_sub
    }
    
    # Figure 5: add a true blank 4th panel
    if (i == 5 && length(top_panels) == 3) {
      blank_top <- ggplot() +
        theme_void() +
        theme(
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
        )
      
      blank_bottom <- ggplot() +
        theme_void() +
        theme(
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.margin = margin(t = 0, r = 2, b = 2, l = 2)
        )
      
      top_panels[["blank"]] <- blank_top
      bottom_panels[["blank"]] <- blank_bottom
    }
    
    while (length(top_panels) < 4) {
      idx_blank <- paste0("blank_", length(top_panels) + 1)
      top_panels[[idx_blank]] <- ggplot() + theme_void()
      bottom_panels[[idx_blank]] <- ggplot() + theme_void()
    }
    
    top_row <- wrap_plots(top_panels, ncol = 4)
    bottom_row <- wrap_plots(bottom_panels, ncol = 4)
    
    combined_plot[[i]] <- top_row / bottom_row +
      plot_layout(heights = c(7.5, 2.5)) &
      theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 5))
    
    print(combined_plot[[i]])
  }
}

# Save as PDF
pdf_files <- paste0("output/sub_figure/att28_death_KM_0", 1:6, "_age_group.pdf")

for (i in seq_along(combined_plot)) {
  if (i <= 3) {
    w <- 16
    h <- 11
  } else {
    w <- 21
    h <- 4.5
  }
  
  CairoPDF(file = pdf_files[i], width = w, height = h)
  print(combined_plot[[i]])
  dev.off()
}

