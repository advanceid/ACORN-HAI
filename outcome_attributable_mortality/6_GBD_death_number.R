# Clear 
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr, 
                 magrittr)
})

# Set working directory
wd <- "./"
setwd(wd)

# Load data
df_raw <- read.csv("data/IHME-GBD.csv")

df <- df_raw %>%
  mutate(across(where(is.numeric),
                ~ format(.x, scientific = FALSE, trim = TRUE))) %>%
  mutate(syndrome = case_when(
    cause_name == "Lower respiratory infections" ~ "LRI",
    cause_name %in% c("Maternal sepsis and other maternal infections",
                      "Neonatal sepsis and other neonatal infections") ~ "BSI",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    val   = as.numeric(val),
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )

df <- df %>%
  mutate(income_level = case_when(
    # --- Low income ---
    location_name %in% c("Democratic People's Republic of Korea") ~ "Low income",
    
    # --- Lower middle income ---
    location_name %in% c("India", "Indonesia", "Philippines", "Viet Nam",
                         "Myanmar", "Pakistan", "Bangladesh", "Nepal",
                         "Bhutan", "Lao People's Democratic Republic",
                         "Cambodia", "Kyrgyzstan", "Uzbekistan",
                         "Tajikistan", "Timor-Leste", "Sri Lanka") ~ "Lower middle income",
    
    # --- Upper middle income ---
    location_name %in% c("China", "Malaysia", "Thailand", "Kazakhstan",
                         "Turkmenistan", "Azerbaijan", "Armenia",
                         "Georgia", "Mongolia", "Maldives") ~ "Upper middle income",
    
    # --- High income ---
    location_name %in% c("Singapore", "Japan", "Republic of Korea",
                         "Brunei Darussalam", "Seychelles", "Mauritius",
                         "Taiwan (Province of China)") ~ "High income",
    
    TRUE ~ "Unclassified"
  ))


# Aggregate by income × syndrome (LRI and BSI separately)
agg_income <- df %>%
  group_by(income_level, syndrome) %>%
  summarise(
    deaths_est  = sum(val,   na.rm = TRUE),
    deaths_low  = sum(lower, na.rm = TRUE),
    deaths_high = sum(upper, na.rm = TRUE),
    .groups = "drop"
  )

# Aggregate by income × all infections (LRI + BSI combined)
agg_income_total <- df %>%
  filter(syndrome %in% c("LRI", "BSI")) %>%
  group_by(income_level) %>%
  summarise(
    deaths_est  = sum(val,   na.rm = TRUE),
    deaths_low  = sum(lower, na.rm = TRUE),
    deaths_high = sum(upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(syndrome = "All infections") %>%
  select(income_level, syndrome, everything())

# Asia overall × syndrome
agg_asia_syndrome <- df %>%
  group_by(syndrome) %>%
  summarise(
    deaths_est  = sum(val,   na.rm = TRUE),
    deaths_low  = sum(lower, na.rm = TRUE),
    deaths_high = sum(upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(income_level = "Overall") %>%
  select(income_level, everything())

# Asia overall × all infections
agg_asia_total <- df %>%
  filter(syndrome %in% c("LRI", "BSI")) %>%
  summarise(
    deaths_est  = sum(val,   na.rm = TRUE),
    deaths_low  = sum(lower, na.rm = TRUE),
    deaths_high = sum(upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(income_level = "Overall", syndrome = "All infections") %>%
  select(income_level, syndrome, everything())

# Combine everything
agg_all <- bind_rows(
  agg_income,
  agg_income_total,
  agg_asia_syndrome,
  agg_asia_total
) %>%
  arrange(factor(income_level,
                 levels = c("Overall","High income","Upper middle income",
                            "Lower middle income","Low income","Unclassified")),
          syndrome) 

# Save
saveRDS(agg_all, "data/GBD_death.RDS")
###