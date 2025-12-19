# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 survminer,
                 gt,
                 purrr,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_all <- readRDS("data/first28_data.RData")
df_used <- split(df_all, df_all$infection_types)

# Load fonts
loadfonts()

# Fit the survreg model
distributions <- c("weibull", "lognormal", "loglogistic", "gaussian", "extreme")

# 
# different distributions
aic_df_list <- list()

for (i in 1:3) {
  df <- df_used[[i]]
  
  aic_values <- numeric(length(distributions))
  
  for (j in seq_along(distributions)) {
    model <- survreg(Surv(time, event) ~ age_new + sex + 
                       country_region + 
                       hpd_admreason + comorbidities_CCI + 
                       severity_score_scale + icu_hd_ap + 
                       aci_car + ent_thir + ent_car, 
                     data = df, dist = distributions[j])
    
    aic_values[j] <- AIC(model)
  }
  
  aic_df_list[[i]] <- data.frame(
    Distribution = distributions,
    AIC = aic_values
  ) %>%
    rename(!!paste0("AIC", i) := AIC)
}


aic_all <- reduce(aic_df_list, left_join, by = "Distribution")
aic_all <- aic_all %>% arrange(AIC1)

# gt table
table <- aic_all %>%
  gt() %>%
  cols_label(
    AIC1 = md("**AIC<sup>[1]</sup>**"),
    AIC2 = md("**AIC<sup>[2]</sup>**"),
    AIC3 = md("**AIC<sup>[3]</sup>**")
  ) %>%
  cols_width(
    1 ~ px(110),
    2:ncol(aic_all) ~ px(80)
  ) %>%
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
      rows = nrow(aic_all)
    )
  ) %>%
  cols_align(
    align = "center",
    columns = 1:ncol(aic_all)
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
  ) %>%
  tab_source_note(
    source_note = html(
      "Note: [1] VAP, [2] Hospital-acquired BSI, [3] Healthcare-associated BSI.<br>Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, AIC = Akaike Information Criterion."
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10)
    ),
    locations = cells_source_notes()
  ) 

print(table)

# Save
gtsave(data = table, filename = "output/table/distributions_AIC_first28_death.html")
