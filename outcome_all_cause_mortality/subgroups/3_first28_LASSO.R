# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 glmnet,
                 car,
                 gt,
                 purrr,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df_all <- readRDS("data/first28_data.RData")
df_used <- split(df_all, df_all$infection_types)

# LASSO
# For variable selection when there is multicollinearity among the predictors
table_vif <- list()

system.time({
  for (i in 1:3) {
    df <- df_used[[i]]
    
    # For num
    x.factors <- model.matrix(~ df$sex + 
                                df$country_region + df$country_income +
                                df$hpd_admreason + df$icu_hd_ap +
                                df$aci_car + df$ent_thir + df$ent_car)[,-1]
    
    # For factor
    x.factors <- as.matrix(data.frame(x.factors, df$age_new, df$comorbidities_CCI, df$severity_score_scale))
    
    x.factors.df <- data.frame(x.factors)
    dummy_y <- rnorm(nrow(x.factors.df))
    
    # Model
    model <- lm(dummy_y ~ ., data = x.factors.df)
    
    # Variance Inflation Factor (VIF)
    vif_values <- vif(model)
    vif_df <- as.data.frame(vif_values)
    vif_df$vif_values <- round(vif_df$vif_values, 3)
    
    # Creat table
    table_row <- list(
      c("Age", NA),
      c(NA, vif_df$vif_values[13]),
      c("Sex", NA),
      c("Female", "—"),
      c("Male", vif_df$vif_values[1]),
      c("Region", NA),
      c("Eastern Mediterranean Region", "—"),
      c("South-East Asian Region", vif_df$vif_values[2]),
      c("Western Pacific Region", vif_df$vif_values[3]),
      c("World Bank income status", NA),
      c("High income", "—"),
      c("Upper middle income", vif_df$vif_values[4]),
      c("Lower middle income", vif_df$vif_values[5]),
      c("Primary admission reason", NA),
      c("Infectious disease", "—"),
      c("Gastrointestinal disorder", vif_df$vif_values[6]),
      c("Pulmonary disease", vif_df$vif_values[7]),
      c("Others", vif_df$vif_values[8]),
      c("Charlson comorbidity index", NA),
      c(NA, vif_df$vif_values[14]),
      c("Severity score of disease", NA),
      c(NA, vif_df$vif_values[15]),
      c("Admission to ICU/HD at enrollment", NA),
      c("No", "—"),
      c("Yes", vif_df$vif_values[9]),
      c("Carbapenem-resistant Acinetobacter spp.", NA),
      c("Absent", "—"),
      c("Present", vif_df$vif_values[10]),
      c("Third-generation cephalosporin-resistant Enterobacterales", NA),
      c("Absent", "—"),
      c("Present", vif_df$vif_values[11]),
      c("Carbapenem-resistant Enterobacterales", NA),
      c("Absent", "—"),
      c("Present", vif_df$vif_values[12])
    )
    
    # Bind the rows into a matrix
    table_vif[[i]] <- do.call(rbind, table_row)
    
    #
    table_vif[[i]] <- table_vif[[i]] %>%
      as.data.frame() %>%
      mutate(across(1, ~ ifelse(.x == "NA", "", .x)))
    
  }
})

#
table_vif_all <- cbind(table_vif[[1]], 
                       table_vif[[2]][,2], 
                       table_vif[[3]][,2])

colnames(table_vif_all) <- c("Variables", levels(df_all$infection_types))

# gt table
lasso_table <- table_vif_all %>%
  gt() %>%
  sub_missing(
    columns = 1:4,
    missing_text = ""
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels(everything())
  ) %>%
  cols_width(
    1 ~ px(200),
    2 ~ px(60),
    3 ~ px(110),
    4 ~ px(125),
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
    data_row.padding = px(0),
    source_notes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(
      rows = nrow(table_vif_all)
    )
  ) %>%
  cols_align(
    align = "center",
    columns = 2:4
  )  %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1)
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Variance inflation factor</span></sup>**"),
    columns = 2:4
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman", 
      size = px(10),             
      weight = "bold"    
    ),
    locations = cells_column_spanners(spanners = tidyselect::matches("Variance inflation factor")) 
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      columns = 1,
      rows = c(1, 3, 6, 10, 14, 19, 21, 23, 26, 29, 32) 
    )
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
  text_transform(
    locations = cells_body(
      columns = Variables,
      rows = Variables == "Carbapenem-resistant Acinetobacter spp."
    ),
    fn = function(x) html("<b>Carbapenem-resistant <i>Acinetobacter</i> spp.</b>")
  )

lasso_table

# Save
gtsave(lasso_table, filename = "output/table/first28_LASSO_sub.html")
###
