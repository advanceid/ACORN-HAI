# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 survival,
                 gt,
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

# 
cox_zph <- list()

system.time({
  for (i in 1:3) {
    df <- df_used[[i]]
    
    # Test the proportional hazards assumption for select MODEL
    model_cox <- coxph(Surv(time, event) ~ age_new + sex + 
                         country_region + country_income +
                         hpd_admreason + comorbidities_CCI + 
                         severity_score_scale +
                         icu_hd_ap + 
                         aci_car + ent_thir + ent_car, 
                       data = df)
    
    
    # Perform the Cox.zph test
    cox_zph_result <- cox.zph(model_cox)
    
    # Convert the Cox.zph result to a data frame
    cox_zph_df <- as.data.frame(cox_zph_result$table)
    
    cox_zph_df$Variables <- c("Age", "Sex", 
                              "Region", "World Bank income status",
                              "Primary admission reason",
                              "Charlson comorbidity index",
                              "Severity score of disease",
                              "Admission to ICU/HD at enrollment",
                              "Carbapenem-resistant Acinetobacter spp.", 
                              "Third-generation cephalosporin-resistant Enterobacterales", 
                              "Carbapenem-resistant Enterobacterales", 
                              "GLOBAL")
    
    rownames(cox_zph_df) <- NULL
    
    # Reorder columns to match the desired output
    cox_zph_df <- cox_zph_df[,c(4,1:3)]
    
    cox_zph_df$p <- ifelse(cox_zph_df$p < 0.001, "<0.001",
                           ifelse(cox_zph_df$p < 0.01, sprintf("%.3f", cox_zph_df$p), 
                                  sprintf("%.2f", cox_zph_df$p)))
  
    cox_zph[[i]] <- cox_zph_df
    
  }
})

cox_zph_all <- cbind(cox_zph[[1]], 
                     cox_zph[[2]][,-1], 
                     cox_zph[[3]][,-1])

colnames(cox_zph_all) <- c("Variables", 
                           "chisq_vap", "df_vap", "p_vap",
                           "chisq_bsi1", "df_bsi1", "p_bsi1",
                           "chisq_bsi2", "df_bsi2", "p_bsi2")

# gt table
table <- cox_zph_all %>%
  gt() %>%
  fmt_number(columns=c(2, 5, 8),decimals= 3) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels(everything())
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>VAP</span>**"),
    columns = 2:4
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Hospital-acquired BSI</span>**"),
    columns = 5:7
  ) %>%
  tab_spanner(
    label = md("**<span style = 'color:black;'>Healthcare-associated BSI</span>**"),
    columns = 8:10
  ) %>%
  cols_label(
    chisq_vap = md("Chi-sq"),
    df_vap = md("DF"),
    p_vap = md("*p*"),
    chisq_bsi1 = md("Chi-sq"),
    df_bsi1 = md("DF"),
    p_bsi1 = md("*p*"),
    chisq_bsi2 = md("Chi-sq"),
    df_bsi2 = md("DF"),
    p_bsi2 = md("*p*")
  ) %>%
  cols_align(
    align = "center",
    columns = 2:ncol(cox_zph_all)
  ) %>%
  cols_align(
    align = "left",
    columns = 1
  ) %>%
  tab_style(
    style = cell_text(align = "left", v_align = "middle"),
    locations = cells_column_labels(columns = 1)
  ) %>%
  cols_width(
    1 ~ px(250),
    2:ncol(cox_zph_all) ~ px(60)
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
      rows = nrow(cox_zph_all)
    )
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(
      rows = nrow(cox_zph_df)
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10)
    ),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",   
      size = px(10),             
      weight = "bold"           
    ),
    locations = cells_column_spanners(spanners = tidyselect::matches("VAP|BSI"))
  ) %>%
  tab_source_note(
    source_note = "Abbreviations: Chi-sq = Chi-Squared Statistic, DF = Degrees of freedom."
  ) %>%
  tab_style(
    style = cell_text(
      font = "Times New Roman",
      size = px(10)
    ),
    locations = cells_source_notes()
  ) %>%
  text_transform(
    locations = cells_body(
      columns = 1,
      rows = Variables == "Carbapenem-resistant Acinetobacter spp."
    ),
    fn = function(x) html("Carbapenem-resistant <i>Acinetobacter</i> spp.")
  )


print(table)
  
# Save
gtsave(data = table, filename = "output/table/cox_zph_first28_death_sub.html")
