# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 magrittr,
                 survival,
                 survminer,
                 extrafont,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/first28_data.RData")

# Load fonts
loadfonts()

# Adjust time for events at 28 days
df$time[df$time == 28 & df$event == 1] <- 27.999

# Fit an overall survival object
surv_object <- Surv(time = df$time, event = df$event)

# Fit Kaplan-Meier model without stratification
km_fit <- survfit(surv_object ~ 1, data = df)

# Define custom theme settings
custom_theme <- theme_minimal(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(size = 10, color = "black", 
                               family = "Times New Roman"), 
    plot.subtitle = element_text(size = 10, color = "black",
                                 family = "Times New Roman"),
    axis.title.x = element_text(size = 10, color = "black",
                                family = "Times New Roman"),
    axis.title.y = element_text(size = 10, color = "black",
                                family = "Times New Roman", margin = margin(r = 10)),
    plot.title = element_text(size = 10, color = "black",
                              family = "Times New Roman"),
    legend.text = element_text(size = 10, color = "black",
                               family = "Times New Roman")
  )

# Create Kaplan-Meier plot without p-value (no comparison needed)
km_plot <- ggsurvplot(
  km_fit,
  data = df,
  font.family = "Times New Roman",
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "black",
  linetype = 1,
  ggtheme = custom_theme,
  title = " ",
  xlab = "Follow-up time since infection onset (days)",
  ylab = "Survival probability",
  legend.title = " ",
  legend.labs = "All infection syndromes",
  palette = "#ce4f55",
  risk.table.title = "Number at risk",
  tables.theme = theme_cleantable(),
  risk.table.y.text = FALSE,
  font.title = c(10, "plain", "black"),
  font.x = c(10, "plain", "black"),
  font.y = c(10, "plain", "black"),
  risk.table.fontsize = 3.5,
  font.tickslab = c(10, "plain", "black"),
  font.risk.table = c(10, "plain", "black"),
  break.time.by = 7,
  xlim = c(0, 28)
)

# Save figure
cairo_pdf(file = "output/pdf/first28_death_KM_overall.pdf", 
          width = 6, height = 4)
print(km_plot, newpage = FALSE)
dev.off()

