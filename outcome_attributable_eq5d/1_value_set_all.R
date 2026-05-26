# Clear
rm(list = ls())

suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr, tidyr, tibble, gt,
    stringr,
    ggplot2,
    scales,
    Cairo,
    grid,
    patchwork
  )
})

# Working directory & data
wd <- "./"
setwd(wd)

# all_eq: all pathogens combined (for Overall pathogen)
# df_all: list of 6 datasets
all_eq <- readRDS("data/att_eq_all.RData")
df_all <- readRDS("data/att_eq.RData")

# convert display strings only
dot_to_mid_char <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  gsub("(?<=\\d)\\.(?=\\d)", "·", x, perl = TRUE)
}

# Format mean (SD)
format_mean_sd <- function(x) {
  if (length(na.omit(x)) == 0) return("")
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

# Infection types order
infection_levels <- c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")

# Priority pathogen labels (matching df_all order)
res_names <- c("CRA", "3GCRE", "CRE", "CRP", "VRE", "MRSA")

# Country to value-set mapping
value_map <- list(
  list(countries = c("Singapore", "Brunei", "Hong Kong SAR China"),
       text = "Singapore (Singapore/Brunei/Hong Kong SAR China)"),
  list(countries = "Japan", text = "Japan"),
  list(countries = "Taiwan", text = "Taiwan"),
  list(countries = "China", text = "China"),
  list(countries = "Malaysia", text = "Malaysia"),
  list(countries = c("Thailand", "Mongolia", "Indonesia"),
       text = "Thailand (Thailand/Mongolia/Indonesia)"),
  list(countries = "Sri Lanka", text = "Sri Lanka"),
  list(countries = c("Bangladesh", "Nepal", "Cambodia", "Philippines",
                     "Pakistan", "Timor-Leste", "Vietnam", "India"),
       text = "Pakistan (Pakistan/Bangladesh/Nepal/Cambodia/Philippines/Timor-Leste/Vietnam/India)")
)

# Collect value-set strings for a set of countries (unique, preserve order)
collect_value_sets <- function(cntrs) {
  texts <- vapply(value_map, function(block) {
    if (any(block$countries %in% cntrs)) block$text else NA_character_
  }, FUN.VALUE = character(1))
  texts <- texts[!is.na(texts)]
  texts <- unique(texts)
  if (length(texts) == 0) "" else paste(texts, collapse = ", ")
}

# Income groups (used for blocks)
income_groups <- if (is.factor(df_all[[1]]$country_income)) {
  levels(df_all[[1]]$country_income)
} else {
  unique(df_all[[1]]$country_income)
}

# All countries present in dataset
all_countries <- Reduce(union, lapply(df_all, function(x) unique(x$country)))

# Block builder
make_block <- function(label, filter_fun, value_text) {
  rows <- c("Overall", infection_levels)
  row_df <- tibble(`Infection type` = rows)
  
  # Priority pathogens
  for (i in seq_along(df_all)) {
    df <- filter_fun(df_all[[i]])
    
    val_overall <- format_mean_sd(df$eq_5d_3l)
    
    vals_by_inf <- df %>%
      group_by(infection_types) %>%
      summarise(val = format_mean_sd(eq_5d_3l), .groups = "drop") %>%
      complete(infection_types = infection_levels, fill = list(val = "")) %>%
      arrange(factor(infection_types, levels = infection_levels)) %>%
      pull(val)
    
    row_df[[res_names[i]]] <- c(val_overall, vals_by_inf)
  }
  
  # Overall pathogen (use all_eq)
  df_all_overall <- filter_fun(all_eq)
  val_overall_all <- format_mean_sd(df_all_overall$eq_5d_3l)
  
  vals_by_inf_all <- df_all_overall %>%
    group_by(infection_types) %>%
    summarise(val = format_mean_sd(eq_5d_3l), .groups = "drop") %>%
    complete(infection_types = infection_levels, fill = list(val = "")) %>%
    arrange(factor(infection_types, levels = infection_levels)) %>%
    pull(val)
  
  row_df[["Overall pathogen"]] <- c(val_overall_all, vals_by_inf_all)
  
  # Add value set description (only on the first row of each block)
  row_df$`Value set used` <- c(value_text, rep("", nrow(row_df) - 1))
  
  # Section header row
  section_row <- tibble(
    `Infection type` = label,
    `Overall pathogen` = "",
    CRA = "",
    `3GCRE` = "",
    CRE = "",
    CRP = "",
    VRE = "",
    MRSA = "",
    `Value set used` = ""
  )
  
  bind_rows(section_row, row_df)
}

# Build summary table
summary_rows <- list()

# 1) All countries/regions
all_value_text <- collect_value_sets(all_countries)
summary_rows[["All countries/regions"]] <- make_block(
  label = "All countries/regions",
  filter_fun = function(df) df,
  value_text = all_value_text
)

# 2) Income-specific blocks
for (income in income_groups) {
  filter_fun <- function(df) dplyr::filter(df, country_income == income)
  
  countries_in_block <- Reduce(union, lapply(df_all, function(x) {
    unique(dplyr::filter(x, country_income == income)$country)
  }))
  
  value_text <- collect_value_sets(countries_in_block)
  
  summary_rows[[income]] <- make_block(
    label = income,
    filter_fun = filter_fun,
    value_text = value_text
  )
}

summary_tbl <- bind_rows(summary_rows) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))

# replace decimal points in display strings only
summary_tbl <- summary_tbl %>%
  mutate(across(everything(), dot_to_mid_char))

# Ensure column order (spanner needs contiguous columns)
summary_tbl <- summary_tbl %>%
  dplyr::select(
    `Infection type`, `Overall pathogen`, CRA, `3GCRE`, CRE, CRP, VRE, MRSA, `Value set used`
  )

# gt table
section_labels <- c("All countries/regions", income_groups)

gt_table <- summary_tbl %>%
  gt() %>%
  tab_spanner(
    label = "Priority pathogens",
    columns = c(CRA, `3GCRE`, CRE, CRP, VRE, MRSA),
    id = "mu"
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
    footnotes.padding = px(0)
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
    locations = cells_body(rows = nrow(summary_tbl))
  ) %>%
  cols_width(
    1 ~ px(130),
    c(2:3, 5:7) ~ px(90),
    c(4, 8) ~ px(110),
    9 ~ px(420)
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(rows = `Infection type` %in% section_labels)
  ) %>%
  cols_align(align = "center", columns = 2:ncol(summary_tbl)) %>%
  cols_align(align = "left", columns = 1) %>%
  tab_style(
    style = cell_text(font = c("Times New Roman"), size = px(10)),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body()
    )
  ) %>%
  tab_style(
    style = cell_text(font = c("Times New Roman"), size = px(10)),
    locations = cells_footnotes()
  ) %>%
  tab_style(
    style = cell_text(font = c("Times New Roman"), size = px(10)),
    locations = cells_source_notes()
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman", size = px(10), weight = "bold"),
    locations = cells_column_spanners(spanners = "mu")
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold", v_align = "middle")),
    locations = cells_column_labels(columns = everything())
  ) %>%
  cols_label(
    `Infection type` = md("Infection syndromes"),
    `Overall pathogen` = md("Overall pathogens<sup>#</sup>"),
    CRA = md("CRA & CSA<sup>#</sup>"),
    `3GCRE` = md("3GCRE & 3GCSE<sup>#</sup>"),
    CRE = md("CRE & CSE<sup>#</sup>"),
    CRP = md("CRP & CSP<sup>#</sup>"),
    VRE = md("VRE & VSE<sup>#</sup>"),
    MRSA = md("MRSA & MSSA<sup>#</sup>"),
    `Value set used` = md("Value set used")
  ) %>%
  tab_source_note(
    html(
      "<sup>#</sup> Values were EQ-5D-3L mean utility (standard deviation).<br>
      Abbreviations: VAP = ventilator-associated pneumonia; BSI = bloodstream infection; 
      CRA = carbapenem-resistant <i>Acinetobacter</i> spp.; CSA = carbapenem-susceptible <i>Acinetobacter</i> spp.; 
      3GCRE = third-generation cephalosporin-resistant Enterobacterales; 3GCSE = third-generation cephalosporin-susceptible Enterobacterales; 
      CRE = carbapenem-resistant Enterobacterales; CSE = carbapenem-susceptible Enterobacterales; 
      CRP = carbapenem-resistant <i>Pseudomonas</i> spp.; CSP = carbapenem-susceptible <i>Pseudomonas</i> spp.; 
      VRE = vancomycin-resistant <i>Enterococcus</i> spp.; VSE = vancomycin-susceptible <i>Enterococcus</i> spp.; 
      MRSA = methicillin-resistant <i>Staphylococcus aureus</i>; MSSA = methicillin-susceptible <i>Staphylococcus aureus</i>."
    )
  )

# Save
gtsave(gt_table, filename = "output/value_set.html")
###
