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

#
df$time[df$time == 28 & df$event == 1] <- 27.999
# Fit a survival object
surv_object <- Surv(time = df$time, event = df$event)

# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ infection_types, data = df)

# p
p_model <- survdiff(surv_object ~ infection_types, data = df)
p_value <- 1 - pchisq(p_model$chisq, df = length(p_model$n) - 1)
#
p_label <- ifelse(p_value < 0.001, "< 0.001",
                  ifelse(p_value < 0.01, paste0("= ", sprintf("%.3f", p_value)),
                         paste0("= ", sprintf("%.2f", p_value))))

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
                                family = "Times New Roman"),
    plot.title = element_text(size = 10, color = "black",
                              family = "Times New Roman"),
    legend.title = element_text(size = 10, color = "black",
                                family = "Times New Roman"),
    legend.text = element_text(size = 10, color = "black",
                               family = "Times New Roman")
  )

# Create Kaplan-Meier plot with custom theme
km_plot <- ggsurvplot(
  km_fit,
  data = df,
  pval = paste0("p ", p_label),
  pval.size = 3.5,
  font.family = "Times New Roman",
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = 1,
  ggtheme = custom_theme,
  title = " ",
  xlab = "Follow-up time since infection onset (days)",
  ylab = "Survival probability",
  legend.title = "",
  legend.labs = c("VAP", 
                  "Hospital-acquired BSI", 
                  "Healthcare-associated BSI"),
  palette = c("#1A75BB","#009344","#F05A28"),
  risk.table.title = "Number at risk",
  tables.theme = theme_cleantable(),
  risk.table.y.text = F,
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
cairo_pdf(file = "output/pdf/first28_death_KM.pdf", 
         width = 6, height = 4)
print(km_plot, newpage = FALSE)
dev.off()
###
