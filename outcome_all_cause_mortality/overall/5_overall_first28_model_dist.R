# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 survminer,
                 gt,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/first28_data.RData")

# Load fonts
loadfonts()

# Fit the survreg model
distributions <- c("weibull", "lognormal", "loglogistic", "gaussian", "extreme")

# 
aic_values <- numeric(length(distributions))

# different distributions
for (i in 1:length(distributions)) {
  model <- survreg(Surv(time, event) ~ age_new + sex + 
                     country_region + country_income +
                     hpd_admreason + comorbidities_CCI + 
                     severity_score_scale + icu_hd_ap + 
                     aci_car + ent_thir + ent_car +
                     infection_types, 
                   data = df, dist = distributions[i])
  aic_values[i] <- AIC(model)
}

# Best model
best_dist <- distributions[which.min(aic_values)]
cat("Best model:", best_dist, "AIC:", min(aic_values), "\n")

# Create a data frame
aic_df <- data.frame(
  Distribution = distributions,
  AIC = aic_values
)

# Sort the data frame by AIC values
aic_df_sorted <- aic_df[order(aic_df$AIC), ]

rownames(aic_df_sorted) <- NULL

table <- aic_df_sorted %>%
  gt() %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    column_labels.border.top.color = "black",
    column_labels.border.top.width = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2),
    table_body.hlines.color = "white",
    table_body.hlines.width = px(0),
    table.border.bottom.color = "white",
    table.border.bottom.width = px(0),
    data_row.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(
      rows = nrow(aic_df_sorted)
    )
  ) %>%
  cols_align(
    align = "center",
    columns = 1:2
  ) %>%
  tab_style(
    style = cell_text(
      font = c("Times New Roman"),
      size = px(10)
    ),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  ) 

#
print(table)

# Save
gtsave(data = table, filename = "output/table/distributions_AIC_first28_death_all.html")
