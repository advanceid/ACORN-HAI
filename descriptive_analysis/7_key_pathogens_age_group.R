# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr,
    dplyr,
    tidyr,
    stringr,
    purrr,
    ggplot2,
    forcats,
    Cairo,
    extrafont,
    tibble,
    grid,
    cowplot,
    ggtext,
    scales
  )
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df  <- readRDS("data/clean_data/baseline_outcomes_index.RData")
ast <- readRDS("data/clean_data/ast_all_index.RData")

# Recode age group
age_levels_plot <- c(
  "<1 year",
  "1–4 years",
  "5–14 years",
  "15–49 years",
  "50–69 years",
  "≥70 years"
)

df$age_group_new <- factor(
  df$age_group_new,
  levels = 1:6,
  labels = age_levels_plot
)

# Add org_combined from AST
org_com_unique <- ast %>%
  select(recordid, org_combined) %>%
  distinct(recordid, .keep_all = TRUE)

df_used <- df %>%
  left_join(org_com_unique, by = "recordid")


# Episode-organism long table
df_pathogen <- df_used %>%
  select(recordid, infection_types, age_group_new, first28_death, org_combined) %>%
  distinct()

df_pathogen_long <- df_pathogen %>%
  filter(!is.na(org_combined), !is.na(age_group_new)) %>%
  mutate(org_combined = str_split(org_combined, pattern = ",\\s*")) %>%
  unnest(org_combined) %>%
  mutate(org_combined = str_trim(org_combined)) %>%
  filter(org_combined != "")


# Join AST resistance variables
df_pathogen_long2 <- df_pathogen_long %>%
  mutate(org_key = str_to_lower(str_squish(org_combined)))

ast_long2 <- ast %>%
  transmute(
    recordid,
    org_key = str_to_lower(str_squish(org_names_all)),
    ris_Carbapenem,
    `ris_Third-generation cephalosporin`,
    Methicillin,
    Vancomycin,
    ris_Fluoroquinolone
  ) %>%
  distinct()

df_pathogen_long_res <- df_pathogen_long2 %>%
  left_join(ast_long2, by = c("recordid", "org_key")) %>%
  select(
    recordid, infection_types, age_group_new, first28_death, org_combined,
    ris_Carbapenem, `ris_Third-generation cephalosporin`,
    Methicillin, Vancomycin, ris_Fluoroquinolone
  )


# Define Enterobacterales
ent_spec_org <- c(
  "E. coli", "K. pneumoniae", "Klebsiella",
  "Enterobacter", "Serratia", "Proteus", "Morganella",
  "Citrobacter", "Providencia", "Raoultella", "Kluyvera",
  "Cronobacter", "Leclercia", "Pluralibacter",
  "Edwardsiella", "Pantoea", "Escherichia", "Salmonella",
  "Kalamiella"
)


##########################
# Fluoroquinolone, Stenotrophomonas spp.
# Burkholderia CR/CS
# VRE, VSE
# MRSA, MSSA
# CRA, CSA
# 3GCSE
# Candida spp.
##########################
# Create mutually exclusive pathogen groups
df_plot <- df_pathogen_long_res %>%
  mutate(
    org_clean = str_squish(org_combined),
    Enterobacterales = ifelse(org_clean %in% ent_spec_org, 1, 0),
    
    pathogen_group = case_when(
      org_clean == "Acinetobacter spp." & ris_Carbapenem == 1 ~ "CRA",
      Enterobacterales == 1 & ris_Carbapenem == 1 ~ "CRE",
      Enterobacterales == 1 & `ris_Third-generation cephalosporin` == 1 ~ "3GCRE",
      org_clean == "Pseudomonas spp." & ris_Carbapenem == 1 ~ "CRP",
      org_clean == "Enterococcus spp." & Vancomycin == 1 ~ "VRE",
      org_clean == "S. aureus" & Methicillin == 1 ~ "MRSA",
      
      org_clean == "Acinetobacter spp." & ris_Carbapenem != 1 ~ "CSA",
      Enterobacterales == 1 & `ris_Third-generation cephalosporin` != 1 ~ "3GCSE",
      org_clean == "Pseudomonas spp." & ris_Carbapenem != 1 ~ "CSP",
      org_clean == "Enterococcus spp." & Vancomycin != 1 ~ "VSE",
      org_clean == "S. aureus" & Methicillin != 1 ~ "MSSA",
      
      org_clean == "Burkholderia" & ris_Carbapenem == 1 ~ "CRB",
      org_clean == "Burkholderia" & ris_Carbapenem != 1 ~ "CSB",
      
      org_clean == "Stenotrophomonas spp." & ris_Fluoroquinolone == 1 ~ "FQ-R Steno. spp.",
      org_clean == "Stenotrophomonas spp." & ris_Fluoroquinolone != 1 ~ "FQ-S Steno. spp.",
      
      org_clean == "Candida spp." ~ "Candida spp.",
      
      TRUE ~ "Others"
    )
  ) %>%
  mutate(
    age_group_new = factor(age_group_new, levels = age_levels_plot),
    pathogen_group = factor(
      pathogen_group,
      levels = c("CRA", "CSA", "CRE",
                 "3GCRE", "3GCSE",
                 "CRP", "CSP",
                 "VRE", "VSE", 
                 "MRSA", "MSSA",
                 "CRB", "CSB", 
                 "FQ-R Steno. spp.", "FQ-S Steno. spp.",
                 "Candida spp.",
                 "Others")
    )
  )


# Count and calculate proportion within each age group
plot_df <- df_plot %>%
  count(age_group_new, pathogen_group, name = "n") %>%
  complete(
    age_group_new = factor(age_levels_plot, levels = age_levels_plot),
    pathogen_group = factor(
      c("CRE","3GCRE", "3GCSE",
        "CRA", "CSA",
        "CRP", "CSP",
        "VRE", "VSE", 
        "MRSA", "MSSA",
        "CRB", "CSB", 
        "FQ-R Steno. spp.", "FQ-S Steno. spp.",
        "Candida spp.",
        "Others"),
      levels = c("CRE","3GCRE", "3GCSE",
                 "CRA", "CSA",
                 "CRP", "CSP",
                 "VRE", "VSE", 
                 "MRSA", "MSSA",
                 "CRB", "CSB", 
                 "FQ-R Steno. spp.", "FQ-S Steno. spp.",
                 "Candida spp.",
                 "Others")
    ),
    fill = list(n = 0)
  ) %>%
  group_by(age_group_new) %>%
  mutate(
    total = sum(n),
    prop = ifelse(total > 0, n / total, 0),
    label = ifelse(prop >= 0.03, scales::percent(prop, accuracy = 1), NA_character_)
  ) %>%
  ungroup()

plot_df$age_group_new <- forcats::fct_rev(plot_df$age_group_new)

# Check each age group sums to 100%
check_sum <- plot_df %>%
  group_by(age_group_new) %>%
  summarise(sum_prop = sum(prop), .groups = "drop")

print(check_sum)

# Define the specific alignment levels (3 items per column logic)
plot_levels_fixed <- c(
  "CRE", "3GCRE", "3GCSE",                       
  "CRA", "CSA", "space1",                       
  "CRP", "CSP", "space2",                       
  "VRE", "VSE", "space3",                    
  "MRSA", "MSSA", "space4",                     
  "CRB", "CSB", "space5",                      
  "FQ-R Steno. spp.", "FQ-S Steno. spp.", "space6", 
  "Candida spp.", "Others", "space7"      
)

# Apply factors to plot_df
plot_df$pathogen_group <- factor(plot_df$pathogen_group, levels = plot_levels_fixed)
plot_df$age_group_new <- forcats::fct_rev(factor(plot_df$age_group_new, levels = age_levels_plot))


# Color
col_ent   <- "#b22222"
col_aci   <- "#Fcae1e"
col_pse   <- "#F3DE8A"
col_enc   <- "#4f7942"
col_sau   <- "#5f8a8b" 
col_bur   <- "#4683b7" 
col_ste   <- "#9e7bb5" 
col_can   <- "#C4CADA" 
col_oth   <- "#D1D1D1"

pathogen_cols <- c(
  "CRE"   = col_ent, "3GCRE" = alpha(col_ent, 0.85), "3GCSE" = alpha(col_ent, 0.7),
  "CRA"   = col_aci, "CSA"   = alpha(col_aci, 0.85),
  "CRP"   = col_pse, "CSP"   = alpha(col_pse, 0.7),
  "VRE"   = col_enc, "VSE"   = alpha(col_enc, 0.85),
  "MRSA"  = col_sau, "MSSA"  = alpha(col_sau, 0.85),
  "CRB"   = col_bur, "CSB"   = alpha(col_bur, 0.85),
  "FQ-R Steno. spp." = col_ste, "FQ-S Steno. spp." = alpha(col_ste, 0.85),
  "Candida spp."     = col_can, 
  "Others"           = col_oth,
  # Placeholders
  "space1"="transparent", "space2"="transparent", "space3"="transparent", 
  "space4"="transparent", "space5"="transparent", "space6"="transparent", "space7"="transparent"
)

# Create a named vector for labels to handle italics and spaces
pathogen_labels <- c(
  "CRE"="CRE", "3GCRE"="3GCRE", "3GCSE"="3GCSE",
  "CRA"="CRA", "CSA"="CSA",
  "CRP"="CRP", "CSP"="CSP",
  "VRE"="VRE", "VSE"="VSE",
  "MRSA"="MRSA", "MSSA"="MSSA",
  "CRB"="CRB", "CSB"="CSB",
  "FQ-R Steno. spp." = expression("FQ-R"~italic("Steno.")~"spp."),
  "FQ-S Steno. spp." = expression("FQ-S"~italic("Steno.")~"spp."),
  "Candida spp."     = expression(italic("Candida")~"spp."),
  "Others"="Others",
  "space1"="", "space2"="", "space3"="", "space4"="", "space5"="", "space6"="", "space7"=""
)

# Final Plotting
p <- ggplot(plot_df, aes(x = prop, y = age_group_new, fill = pathogen_group)) +
  geom_col(
    width = 0.8,
    color = NA,
    linewidth = 0.3,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    family = "Times New Roman",
    size = 3.2,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = pathogen_cols,
    labels = pathogen_labels,
    drop = FALSE # Crucial: prevents the removal of empty "space" levels
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0),
    oob = scales::squish
  ) +
  labs(x = "Proportion of pathogens", y = NULL, fill = "Pathogens") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.2),
    axis.ticks.x = element_line(color = "black", linewidth = 0.2),
    
    # Legend Customization
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    legend.key.width = unit(1, "lines"),
    legend.key.height = unit(0.8, "lines"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.spacing.y = unit(0, "cm"),

    legend.key = element_rect(fill = NA, color = NA),
    
    plot.margin = margin(t = 5, r = 25, b = 5, l = 10, unit = "pt")
  ) +
  # Use guide_legend to force the grid structure
  guides(fill = guide_legend(
    ncol = 8,         
    byrow = FALSE,    
    title.position = "top",
    override.aes = list(color = "white") 
  ))

print(p)

# Save
CairoPDF(
  file = "output/figure/pathogen_age_proportion.pdf",
  width = 7.5,
  height = 4
)
print(p)
dev.off()
###
