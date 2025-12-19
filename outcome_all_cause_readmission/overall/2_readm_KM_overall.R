# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 tidycmprsk, 
                 ggsurvfit,
                 ggplot2,
                 grid,
                 patchwork,
                 tidyverse,
                 survminer,
                 survival,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/readmission90_data.RData")

# Load fonts
loadfonts()

# Convert outcome variable to numeric
df$readm_event <- as.numeric(df$readm_event)

#
df$readm_death_time[df$readm_death_time == 90 & df$readm_event %in% c(1,2)] <- 89.999

# Set factor levels and labels for readm_event
df$readm_event <- factor(df$readm_event,
                         levels = c(0, 2, 1),
                         labels = c("censor", "All-cause death", "Readmission"))

# Cumulative incidence
cuminc <- cuminc(Surv(readm_death_time, readm_event) ~ 1, data = df)
cuminc_data <- tidy_cuminc(cuminc)

cuminc_data <- cuminc_data %>%
  mutate(n.risk = first(n.risk) - cumsum(n.event)) %>%  
  ungroup()

# Customize gg_cuminc plot
plot_data <- cuminc_data %>%
  filter(outcome %in% c("Readmission", "All-cause death")) %>%
  mutate(outcome = factor(outcome, levels = c("Readmission", "All-cause death")))

gg_cuminc <- ggplot(plot_data, aes(x = time, y = estimate,
                                   color = outcome, linetype = outcome)) +
  geom_step(size = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = outcome),
              alpha = 0.2, color = NA) +
  scale_color_manual(
    name = "Events",
    values = c("Readmission" = "#ce4f55", "All-cause death" = "#ce4f55")
  ) +
  scale_fill_manual(
    name = "Events",
    values = c("Readmission" = "#ce4f55", "All-cause death" = "#ce4f55")
  ) +
  scale_linetype_manual(
    name = "Events",
    values = c("Readmission" = "solid", "All-cause death" = "twodash")
  ) +
  scale_x_continuous(breaks = seq(0, 90, by = 15), limits = c(0, 90)) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, family = "Times New Roman"),
    legend.text = element_text(size = 10, family = "Times New Roman"),
    axis.text = element_text(size = 10, family = "Times New Roman"),
    axis.title = element_text(size = 10, family = "Times New Roman"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 17)
  ) +
  labs(x = "Follow-up time (days)", y = "Cumulative incidence") +
  annotate("text", x = 4.1, y = 0.3, label = "p < 0.001 (all-cause death)",
           size = 3.5, color = "black", family = "Times New Roman") +
  annotate("text", x = 3.5, y = 0.33, label = "p < 0.001 (readmission)",
           size = 3.5, color = "black", family = "Times New Roman")


# Build the risk table plot
gg_risktable <- 
  cuminc_data |>
  filter(time %in% seq(0, 90, by = 15)) |>
  select(outcome, time, n.risk, cum.event) %>%
  {
    mutate(., outcome = "Atrisk", stat = n.risk) |> 
      select(outcome, time, stat) |>
      bind_rows(
        select(., outcome = outcome, time, stat = cum.event)
      )
  } |>
  mutate(outcome = recode(outcome, 
                          "All-cause death" = "Adeath", 
                          "Readmission" = "Readm", 
                          "Atrisk" = "Atrisk")) |> 
  mutate(outcome = factor(outcome, levels = c("Adeath", "Readm", "Atrisk"))) |>
  ggplot(aes(x = time, y = factor(outcome), label = stat, linetype = outcome)) +  
  geom_text(size = 3.5, color = "black", family = "Times New Roman",
            check_overlap = TRUE) +
  labs(y = NULL, x = NULL) +
  scale_x_continuous(breaks = seq(0, 90, by = 15), limits = c(0, 90)) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.text = element_text(size = 0),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", family = "Times New Roman"),
    plot.margin = margin(t = 10, r = 5, b = 5, l = 35)
  ) +
  scale_linetype_manual(
    name = "Events", 
    values = c("Adeath" = "dashed", "Readm" = "solid", "Atrisk" = "solid"),  
    breaks = c("Adeath", "Readm", "Atrisk"),  
    labels = c("All-cause death", "Readmission", "At risk")
  )

# Combine
final_plot <- gg_cuminc / gg_risktable + patchwork::plot_layout(heights = c(3, 1))
print(final_plot)

# Save
cairo_pdf(file = "output/pdf/readmission_CIF_overall.pdf", 
          width = 12, height = 5)
print(final_plot)
dev.off()
###