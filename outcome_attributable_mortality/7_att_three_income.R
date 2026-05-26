# Clear & packages
rm(list = ls())

suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    dplyr, magrittr, tidyr, purrr, tibble, gt,
    stringr, htmltools
  )
})

# Working directory & data
wd <- "./"
setwd(wd)

paf   <- readRDS("data/att28_paf_marginal.RData")
death <- readRDS("data/GBD_death.RDS")

# Settings
pathogens <- c("mdr_gnb", "mdr", "amr")

# Keep needed rows
paf <- paf %>%
  filter(
    pathogen %in% pathogens,
    subgroup %in% c("overall", "country_income")
  )

death <- death %>%
  filter(syndrome == "All infections")

# Helpers ---------------------------------------------------------------

num_2 <- function(x, digits = 2) {
  ifelse(
    is.na(x),
    NA,
    formatC(round(x, digits), format = "f", digits = digits)
  )
}

fmt_ci_number <- function(est, lo, hi, digits = 0) {
  est_f <- ifelse(
    is.na(est),
    NA,
    formatC(round(est, digits), format = "f", digits = digits)
  )
  lo_f <- ifelse(
    is.na(lo),
    NA,
    formatC(round(lo, digits), format = "f", digits = digits)
  )
  hi_f <- ifelse(
    is.na(hi),
    NA,
    formatC(round(hi, digits), format = "f", digits = digits)
  )
  
  ifelse(
    is.na(lo) | is.na(hi),
    est_f,
    paste0(est_f, " [", lo_f, ", ", hi_f, "]")
  )
}

fmt_ci_percent_mix <- function(est, lo, hi, digits = 2) {
  est_f <- num_2(100 * est, digits)
  lo_f  <- num_2(100 * lo, digits)
  hi_f  <- num_2(100 * hi, digits)
  
  ifelse(
    is.na(lo) | is.na(hi),
    est_f,
    paste0(est_f, " [", lo_f, ", ", hi_f, "]")
  )
}

# only convert display strings after all calculations are done
dot_to_mid_char <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  gsub("(?<=\\d)\\.(?=\\d)", "·", x, perl = TRUE)
}

# robust parser: extract first 3 numbers from each PAF_CI string
parse_paf_ci <- function(x) {
  nums <- stringr::str_extract_all(x, "-?\\d+\\.?\\d*")
  
  tibble(
    PAF_overall    = purrr::map_dbl(nums, ~ if (length(.x) >= 1) as.numeric(.x[1]) / 100 else NA_real_),
    PAF_overall_lo = purrr::map_dbl(nums, ~ if (length(.x) >= 2) as.numeric(.x[2]) / 100 else NA_real_),
    PAF_overall_hi = purrr::map_dbl(nums, ~ if (length(.x) >= 3) as.numeric(.x[3]) / 100 else NA_real_)
  )
}

# Parse PAF from file
paf2 <- paf %>%
  bind_cols(parse_paf_ci(.$PAF_CI))

# Convert PAF table to income-level table
paf_by_income <- paf2 %>%
  mutate(
    income_level = case_when(
      subgroup == "overall" & level == "overall" ~ "Overall",
      subgroup == "country_income" ~ level,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(income_level)) %>%
  select(
    income_level, pathogen,
    PAF_overall, PAF_overall_lo, PAF_overall_hi
  ) %>%
  distinct()

print(paf_by_income)

# Normalize death labels
death <- death %>%
  mutate(
    income_level = gsub("-", " ", income_level),
    income_level = trimws(income_level)
  )

# If Overall missing, create it
if (!("Overall" %in% unique(death$income_level))) {
  death_overall <- death %>%
    group_by(syndrome) %>%
    summarise(
      deaths_est  = sum(deaths_est,  na.rm = TRUE),
      deaths_low  = sum(deaths_low,  na.rm = TRUE),
      deaths_high = sum(deaths_high, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(income_level = "Overall") %>%
    relocate(income_level)
  
  death <- bind_rows(death, death_overall)
}

# Keep only needed income levels
death <- death %>%
  filter(income_level %in% c(
    "Overall", "High income", "Upper middle income",
    "Lower middle income", "Low income"
  ))

# Expand death rows to all pathogens
death_expanded <- death %>%
  select(income_level, syndrome, deaths_est, deaths_low, deaths_high) %>%
  tidyr::crossing(tibble(pathogen = pathogens))

# Join with PAF
death_joined <- death_expanded %>%
  left_join(paf_by_income, by = c("income_level", "pathogen"))

# Calculate attributable deaths
death_attr_income <- death_joined %>%
  mutate(
    attrib_deaths_est = deaths_est  * PAF_overall,
    attrib_deaths_lo  = deaths_low  * PAF_overall_lo,
    attrib_deaths_hi  = deaths_high * PAF_overall_hi
  ) %>%
  arrange(
    match(income_level, c("Overall", "High income", "Upper middle income", "Lower middle income", "Low income")),
    match(pathogen, pathogens)
  ) %>%
  filter(income_level != "Low income") %>%
  filter(!is.na(income_level), !is.na(pathogen))

# Clean table output
table_clean <- death_attr_income %>%
  mutate(
    Deaths_raw = fmt_ci_number(deaths_est, deaths_low, deaths_high, digits = 0),
    PAF_overall_fmt_raw = fmt_ci_percent_mix(PAF_overall, PAF_overall_lo, PAF_overall_hi, digits = 2),
    Attrib_deaths_raw = fmt_ci_number(attrib_deaths_est, attrib_deaths_lo, attrib_deaths_hi, digits = 0),
    
    Deaths = dot_to_mid_char(Deaths_raw),
    PAF_overall_fmt = dot_to_mid_char(PAF_overall_fmt_raw),
    Attrib_deaths = dot_to_mid_char(Attrib_deaths_raw),
    
    income_level = factor(
      income_level,
      levels = c("Overall", "High income", "Upper middle income", "Lower middle income")
    ),
    pathogen = factor(pathogen, levels = pathogens)
  ) %>%
  filter(!is.na(income_level), !is.na(pathogen)) %>%
  arrange(income_level, pathogen) %>%
  select(
    income_level, pathogen,
    Deaths,
    PAF_overall = PAF_overall_fmt,
    Attrib_deaths
  )

# For table
tbl0 <- table_clean %>%
  filter(!is.na(income_level)) %>%
  mutate(
    pathogen = factor(
      recode(
        as.character(pathogen),
        "amr"     = "AMR",
        "mdr"     = "MDR",
        "mdr_gnb" = "MDR-GNB"
      ),
      levels = c("AMR", "MDR", "MDR-GNB")
    ),
    income_level = factor(
      income_level,
      levels = c("Overall", "High income", "Upper middle income", "Lower middle income")
    )
  ) %>%
  arrange(income_level, pathogen) %>%
  rename(
    `Pathogens (all-cause deaths)` = pathogen,
    `All-cause deaths` = Deaths,
    `Marginal population attributable fractions (%, 95%CI)` = PAF_overall,
    `Attributable deaths (95%CI)` = Attrib_deaths
  )

# Build group labels
group_labels <- tbl0 %>%
  filter(!is.na(income_level)) %>%
  group_by(income_level) %>%
  summarise(
    deaths_first = first(`All-cause deaths`),
    .groups = "drop"
  ) %>%
  mutate(
    income_level = as.character(income_level),
    label = paste0(income_level, "<br>(", deaths_first, ")"),
    id = paste0("grp_", row_number())
  )

# Remove all-cause deaths from body
tbl_show_body <- tbl0 %>%
  mutate(`All-cause deaths` = "") %>%
  select(-income_level, -`All-cause deaths`)

last_row <- nrow(tbl_show_body)

# Render gt
tbl <- gt(tbl_show_body) %>%
  {
    gt_obj <- .
    gl <- group_labels
    
    for (i in seq_len(nrow(gl))) {
      row_idx <- which(as.character(tbl0$income_level) == gl$income_level[i])
      if (length(row_idx) > 0) {
        gt_obj <- gt_obj %>%
          tab_row_group(
            label = md(gl$label[i]),
            id = gl$id[i],
            rows = row_idx
          )
      }
    }
    gt_obj
  } %>%
  row_group_order(na.omit(group_labels$id)) %>%
  cols_align(align = "left", columns = "Pathogens (all-cause deaths)") %>%
  cols_align(
    align = "center",
    columns = c(
      "Marginal population attributable fractions (%, 95%CI)",
      "Attributable deaths (95%CI)"
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    column_labels.border.top.color    = "black",
    column_labels.border.top.width    = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(2),
    table_body.hlines.color           = "transparent",
    table_body.hlines.width           = px(0),
    row_group.border.top.color        = "transparent",
    row_group.border.top.width        = px(0),
    row_group.border.bottom.color     = "transparent",
    row_group.border.bottom.width     = px(0),
    table.border.top.color            = "transparent",
    table.border.top.width            = px(0),
    table.border.bottom.color         = "transparent",
    table.border.bottom.width         = px(0),
    data_row.padding                  = px(0),
    footnotes.padding                 = px(0)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = c("top", "bottom"),
      color = "transparent",
      weight = px(0)
    ),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "transparent",
      weight = px(0)
    ),
    locations = cells_body(rows = everything(), columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(rows = last_row, columns = everything())
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman", size = px(10)),
    locations = list(
      cells_title(groups = "title"),
      cells_title(groups = "subtitle"),
      cells_column_labels(),
      cells_body(),
      cells_footnotes(),
      cells_source_notes()
    )
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman", size = px(10), weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  cols_width(
    1 ~ px(150),
    2 ~ px(160),
    3 ~ px(130)
  ) %>%
  tab_source_note(
    source_note = "Abbreviations: AMR = Antimicrobial Resistance, MDR = Multidrug Resistance, GNB = Gram-negative Bacteria, CI = Confidence Interval."
  ) %>%
  tab_style(
    style = cell_text(v_align = "middle"),
    locations = list(
      cells_body(columns = everything()),
      cells_column_labels(columns = everything()),
      cells_row_groups()
    )
  ) %>%
  tab_options(
    footnotes.marks = c("#", "*", "†", "‡", "§"),
    footnotes.sep = " ",
    footnotes.padding = px(0)
  ) %>%
  tab_footnote(
    footnote = "Data are from cause-specific mortality estimates of the Global Burden of Disease.",
    locations = cells_column_labels(columns = "Pathogens (all-cause deaths)")
  )

tbl

# Save
gtsave(tbl, "output/attrib_deaths_asia.html")

