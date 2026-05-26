# Clear
rm(list = ls())

# Load required packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr,
    magrittr,
    gt,
    extrafont,
    purrr,
    stringr,
    openxlsx,
    htmltools
  )
})

# Set working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Create output folder if not exists
if (!dir.exists("output")) dir.create("output")

# Function: format all decimals to 1 dp in strings
format_1dp <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  
  m <- gregexpr("\\-?\\d+(?:\\.\\d+)?", x, perl = TRUE)
  
  regmatches(x, m) <- lapply(regmatches(x, m), function(v) {
    out <- sprintf("%.1f", as.numeric(v))
    gsub("(?<=\\d)\\.(?=\\d)", "·", out, perl = TRUE)
  })
  
  x
}

# Function: make table
make_att28_table <- function(input_file, output_prefix) {
  
  # Load data
  df_marg <- readRDS(input_file)
  
  # Keep only required pathogens
  keep_pathogen <- c(
    "aci_car",
    "ent_thir",
    "ent_car",
    "pse_car",
    "entc_van",
    "sa_meth"
  )
  
  # Pathogen labels
  pathogen_labels <- c(
    aci_car  = "Carbapenem-resistant Acinetobacter spp.",
    ent_thir = "Third-generation cephalosporin-resistant Enterobacterales",
    ent_car  = "Carbapenem-resistant Enterobacterales",
    pse_car  = "Carbapenem-resistant Pseudomonas spp.",
    entc_van = "Vancomycin-resistant Enterococcus spp.",
    sa_meth  = "Methicillin-resistant Staphylococcus aureus"
  )
  
  # Subgroup labels
  subgroup_labels <- c(
    overall          = "Overall",
    infection_types  = "Infection syndromes",
    age_group_new    = "Age group",
    country_income   = "World Bank income status"
  )
  
  # Desired subgroup order
  subgroup_order <- c("overall", "infection_types", "age_group_new", "country_income")
  
  # Orders within subgroup
  infection_order <- c(
    "VAP",
    "Hospital-acquired BSI",
    "Healthcare-associated BSI"
  )
  
  age_order <- c(
    "<1 year",
    "1–4 years", "1-4 years",
    "5–14 years", "5-14 years",
    "15–49 years", "15-49 years",
    "50–69 years", "50-69 years",
    "≥70 years", ">=70 years"
  )
  
  income_order <- c(
    "High income",
    "Upper middle income",
    "Lower middle income"
  )
  
  # Prepare data
  df <- df_marg %>%
    dplyr::filter(pathogen %in% keep_pathogen) %>%
    dplyr::select(pathogen, subgroup, level, num, Prevalence, AR_CI, PAF_CI) %>%
    dplyr::rename(
      AR_marg  = AR_CI,
      PAF_marg = PAF_CI
    ) %>%
    mutate(
      pathogen = recode(pathogen, !!!pathogen_labels),
      pathogen = factor(
        pathogen,
        levels = unname(pathogen_labels[keep_pathogen])
      ),
      subgroup_display = recode(subgroup, !!!subgroup_labels),
      subgroup = factor(subgroup, levels = subgroup_order),
      level = as.character(level),
      level = ifelse(level == "overall", "Overall", level),
      level_ord = case_when(
        subgroup == "overall" ~ 1,
        subgroup == "infection_types" ~ match(level, infection_order),
        subgroup == "age_group_new" ~ match(level, age_order),
        subgroup == "country_income" ~ match(level, income_order),
        TRUE ~ NA_real_
      ),
      level_ord = ifelse(is.na(level_ord), 999, level_ord)
    ) %>%
    arrange(pathogen, subgroup, level_ord, level) %>%
    mutate(pathogen = as.character(pathogen))
  
  # Build display table
  pathogen_vec <- unique(df$pathogen)
  result_list <- list()
  
  for (p in pathogen_vec) {
    
    tmp_p <- df %>% filter(pathogen == p)
    
    # Pathogen title row
    result_list[[length(result_list) + 1]] <- tibble(
      `Key antibiotic resistant pathogens` = p,
      ` ` = "",
      `Proportion of cases exposed (%)` = "",
      `Attributable risk (%, 95%CI)` = "",
      `Population attributable fraction (%, 95%CI)` = ""
    )
    
    # Data rows
    tmp_disp <- tmp_p %>%
      mutate(
        left_label = case_when(
          subgroup == "overall" ~ "Overall",
          subgroup == "infection_types" ~ "Infection syndromes",
          subgroup == "age_group_new" ~ "Age group",
          subgroup == "country_income" ~ "World Bank income status",
          TRUE ~ subgroup_display
        ),
        right_label = case_when(
          subgroup == "overall" ~ num,
          TRUE ~ paste0(level, " (", gsub("^N\\s*=\\s*", "", num), ")")
        )
      ) %>%
      transmute(
        `Key antibiotic resistant pathogens` = left_label,
        ` ` = right_label,
        `Proportion of cases exposed (%)` = Prevalence,
        `Attributable risk (%, 95%CI)` = AR_marg,
        `Population attributable fraction (%, 95%CI)` = PAF_marg
      )
    
    result_list[[length(result_list) + 1]] <- tmp_disp
  }
  
  result_all <- bind_rows(result_list)
  
  # Format numbers
  cols_to_format <- c(
    "Proportion of cases exposed (%)",
    "Attributable risk (%, 95%CI)",
    "Population attributable fraction (%, 95%CI)"
  )
  
  result_all[cols_to_format] <- lapply(result_all[cols_to_format], format_1dp)
  
  # Identify pathogen title rows
  pathogen_rows <- which(
    result_all$` ` == "" &
      result_all$`Proportion of cases exposed (%)` == ""
  )
  
  # Save Excel
  write.xlsx(
    result_all,
    file = paste0("output/", output_prefix, ".xlsx"),
    overwrite = TRUE
  )
  
  # Build gt table
  table <- result_all %>%
    gt() %>%
    cols_label(
      ` ` = "",
      `Attributable risk (%, 95%CI)` = md("Attributable risk<br>(%, 95%CI)"),
      `Population attributable fraction (%, 95%CI)` = md("Population attributable fraction<br>(%, 95%CI)")
    ) %>%
    fmt_markdown(columns = c(`Key antibiotic resistant pathogens`)) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_column_labels(everything())
    ) %>%
    tab_style(
      style = list(cell_text(weight = "bold")),
      locations = cells_body(
        rows = pathogen_rows,
        columns = `Key antibiotic resistant pathogens`
      )
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
      footnotes.padding = px(0),
      source_notes.padding = px(0)
    ) %>%
    tab_style(
      style = cell_borders(
        sides = "bottom",
        color = "black",
        weight = px(2)
      ),
      locations = cells_body(
        rows = nrow(result_all)
      )
    ) %>%
    cols_align(
      align = "left",
      columns = c(`Key antibiotic resistant pathogens`, ` `)
    ) %>%
    cols_align(
      align = "center",
      columns = c(
        `Proportion of cases exposed (%)`,
        `Attributable risk (%, 95%CI)`,
        `Population attributable fraction (%, 95%CI)`
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
    tab_style(
      style = cell_text(
        font = c("Times New Roman"),
        size = px(10)
      ),
      locations = cells_footnotes()
    ) %>%
    tab_style(
      style = cell_text(
        font = "Times New Roman",
        size = px(10)
      ),
      locations = cells_source_notes()
    ) %>%
    tab_style(
      style = cell_text(
        weight = "bold",
        font = "Times New Roman",
        size = px(10)
      ),
      locations = cells_column_spanners()
    ) %>%
    cols_width(
      `Key antibiotic resistant pathogens` ~ px(260),
      ` ` ~ px(155),
      `Proportion of cases exposed (%)` ~ px(120),
      `Attributable risk (%, 95%CI)` ~ px(145),
      `Population attributable fraction (%, 95%CI)` ~ px(165)
    ) %>%
    tab_footnote(
      footnote = "The proportions of cases exposed were from the ACORN-HAI study.",
      locations = cells_column_labels(columns = "Proportion of cases exposed (%)")
    ) %>%
    tab_footnote(
      footnote = "Only subgroups with total n ≥ 30 were shown.",
      locations = cells_body(
        columns = "Key antibiotic resistant pathogens",
        rows = grepl(
          "Infection syndromes|Age group|World Bank income status",
          result_all$`Key antibiotic resistant pathogens`
        )
      )
    ) %>%
    tab_source_note(
      source_note = "Abbreviations: VAP = Ventilator-Associated Pneumonia, BSI = Bloodstream Infection, CI = Confidence Interval."
    ) %>%
    tab_style(
      style = cell_text(align = "center", v_align = "middle"),
      locations = cells_column_labels(columns = everything())
    ) %>%
    tab_options(footnotes.marks = c("#", "‡"))
  
  # Italicize organism names in pathogen title rows
  table <- table %>%
    text_transform(
      locations = cells_body(
        columns = `Key antibiotic resistant pathogens`,
        rows = result_all$`Key antibiotic resistant pathogens` == "Carbapenem-resistant Acinetobacter spp."
      ),
      fn = function(x) html("<b>Carbapenem-resistant <i>Acinetobacter</i> spp.</b>")
    ) %>%
    text_transform(
      locations = cells_body(
        columns = `Key antibiotic resistant pathogens`,
        rows = result_all$`Key antibiotic resistant pathogens` == "Carbapenem-resistant Pseudomonas spp."
      ),
      fn = function(x) html("<b>Carbapenem-resistant <i>Pseudomonas</i> spp.</b>")
    ) %>%
    text_transform(
      locations = cells_body(
        columns = `Key antibiotic resistant pathogens`,
        rows = result_all$`Key antibiotic resistant pathogens` == "Vancomycin-resistant Enterococcus spp."
      ),
      fn = function(x) html("<b>Vancomycin-resistant <i>Enterococcus</i> spp.</b>")
    ) %>%
    text_transform(
      locations = cells_body(
        columns = `Key antibiotic resistant pathogens`,
        rows = result_all$`Key antibiotic resistant pathogens` == "Methicillin-resistant Staphylococcus aureus"
      ),
      fn = function(x) html("<b>Methicillin-resistant <i>Staphylococcus aureus</i></b>")
    )
  
  print(table)
  
  # Save HTML
  gtsave(
    data = table,
    filename = paste0("output/", output_prefix, ".html")
  )
  
  invisible(result_all)
}

# Original marginal result
result_original <- make_att28_table(
  input_file = "data/att28_paf_marginal.RData",
  output_prefix = "att28_mar_paf"
)

# Mixed-effects sensitivity analysis
result_mix <- make_att28_table(
  input_file = "data/att28_paf_marginal_mix.RData",
  output_prefix = "att28_mar_paf_mix"
)

# Monomicrobial sensitivity analysis
result_mono <- make_att28_table(
  input_file = "data/att28_paf_marginal_mono.RData",
  output_prefix = "att28_mar_paf_mono"
)


