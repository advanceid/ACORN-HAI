# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr,
                 dplyr,
                 ggplot2,
                 tidyr,
                 extrafont,
                 forcats,
                 Cairo)
})

# Import fonts
loadfonts() 

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df  <- readRDS("data/clean_data_RData/data_table_index_new_delete.RData")
ast <- readRDS("data/clean_data_RData/ast_all_index.RData")

# ENT species
ent_spec_org <- c(
  "E. coli", "K. pneumoniae", "Klebsiella", 
  "Enterobacter", "Serratia", "Proteus", "Morganella",
  "Citrobacter", "Providencia", "Raoultella", "Kluyvera",
  "Cronobacter", "Leclercia", "Pluralibacter", 
  "Edwardsiella", "Pantoea", "Escherichia", "Salmonella",
  "Kalamiella"
)

# Identify pathogen IDs
id_aci  <- unique(ast$recordid[ast$org_names_all == "Acinetobacter spp."])
id_ent  <- unique(ast$recordid[ast$org_names_all %in% ent_spec_org])
id_pse  <- unique(ast$recordid[ast$org_names_all == "Pseudomonas spp."])
id_entc <- unique(ast$recordid[ast$org_names_all == "Enterococcus spp."])
id_sa   <- unique(ast$recordid[ast$org_names_all == "S. aureus"])

# Filter datasets
df_aci  <- df %>% filter(recordid %in% id_aci)
df_ent  <- df %>% filter(recordid %in% id_ent)
df_pse  <- df %>% filter(recordid %in% id_pse)
df_entc <- df %>% filter(recordid %in% id_entc)
df_sa   <- df %>% filter(recordid %in% id_sa)

# Define pathogen list
pathogen_list <- list(
  list(data = df_aci,  var = "aci_car",    name = "CRA"),
  list(data = df_ent,  var = "ent_thir",   name = "3GCRE"),
  list(data = df_ent,  var = "ent_car",    name = "CRE"),
  list(data = df_pse,  var = "pse_car",    name = "CRP"),
  list(data = df_entc, var = "entc_van",   name = "VRE"),
  list(data = df_sa,   var = "sa_meth",    name = "MRSA")
)

# Count deaths function
get_death_counts <- function(df, resistant_var, pathogen) {
  df_split <- split(df, df$infection_types)
  do.call(rbind, lapply(names(df_split), function(type) {
    data_i <- df_split[[type]]
    data.frame(
      pathogens         = pathogen,
      infection_type    = type,
      total_patients    = nrow(data_i),
      total_deaths      = sum(data_i$first28_death == 1, na.rm = TRUE),
      resistant_deaths  = sum(
        grepl("resistant$", data_i[[resistant_var]]) &
          data_i$first28_death == 1,
        na.rm = TRUE
      )
    )
  }))
}

# Apply across all pathogens
final_results <- do.call(rbind, lapply(pathogen_list, function(p) {
  get_death_counts(p$data, p$var, p$name)
})) 

# --------
# Combine CRE + 3GCRE 
cre3gcre_raw <- final_results %>%
  filter(pathogens %in% c("CRE", "3GCRE")) %>%
  mutate(non_resistant_deaths = total_deaths - resistant_deaths) %>%
  select(pathogens, infection_type, total_patients, total_deaths, resistant_deaths)

# total_patients
gray_bar <- cre3gcre_raw %>%
  filter(pathogens == "CRE") %>%
  mutate(pathogens = "Enterobacterales",
         count     = total_patients,
         layer     = "total_patients")

# total_deaths
orange_bar <- cre3gcre_raw %>%
  filter(pathogens == "CRE") %>%
  mutate(pathogens = "Enterobacterales",
         count     = total_deaths,
         layer     = "total_deaths")

# resistant_deaths (CRE vs 3GCRE)
resistant_stack <- cre3gcre_raw %>%
  mutate(layer = case_when(
    pathogens == "CRE"   ~ "res_CRE",
    pathogens == "3GCRE" ~ "res_3GCRE"
  ),
  pathogens = "Enterobacterales",
  count     = resistant_deaths
  ) %>%
  select(pathogens, infection_type, count, layer)

# Combine and force factor order on `layer`
plot_data <- bind_rows(
  gray_bar %>% select(pathogens, infection_type, count, layer),
  orange_bar %>% select(pathogens, infection_type, count, layer),
  resistant_stack
) %>%
  mutate(
    infection_type = forcats::fct_rev(factor(
      infection_type,
      levels = c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")
    )),
    layer = factor(layer,
                   levels = c(
                     "total_patients",
                     "total_deaths",
                     "res_3GCRE",
                     "res_CRE"
                   )
    )
  )

# Plot
p2 <- ggplot() +
  geom_col(
    data = filter(plot_data, layer == "total_patients"),
    aes(y = infection_type, x = count),
    fill = "grey90", width = 0.6, color = NA
  ) +
  geom_col(
    data = filter(plot_data, layer == "total_deaths"),
    aes(y = infection_type, x = count),
    fill = "grey65", width = 0.4, color = NA
  ) +
  geom_col(
    data = filter(plot_data, layer %in% c("res_CRE", "res_3GCRE")),
    aes(y = infection_type, x = count, fill = layer),
    width = 0.2, position = "stack", color = NA
  ) +
  scale_fill_manual(
    values = c(
      "res_3GCRE" = "#cbdde6",
      "res_CRE"   = "#194a55"
    ),
    labels = c("3GCRE", "CRE"),
    name   = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 2500),
    breaks = seq(0, 2500, 500)
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    panel.grid.major.y    = element_blank(),
    plot.title            = element_text(
      hjust = 0.5, size = 12, face = "bold",
      family = "Times New Roman"
    ),
    axis.text.y = element_text(
      size = 12, family = "Times New Roman",
      hjust = 1, margin = margin(r = -10)
    ),
    axis.text.x = element_text(
      size = 12, family = "Times New Roman"
    ),
    legend.position = c(0.8, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(
      fill = alpha("white", 0.8), color = NA
    ),
    legend.box.background = element_rect(color = NA),
    legend.text = element_text(
      size = 12, family = "Times New Roman"
    )
  ) +
  labs(
    title = "Enterobacterales",
    x = NULL,
    y = NULL
  )


print(p2)


# -----
# Plot others
draw_pathogen_plot <- function(pathogen_name, organism_label) {
  plot_data <- final_results %>%
    filter(pathogens == pathogen_name) %>%
    mutate(
      infection_type = forcats::fct_rev(factor(
        infection_type,
        levels = c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")
      ))
    )
  
  gray_bar <- plot_data %>%
    mutate(layer = "total_patients", count = total_patients)
  
  orange_bar <- plot_data %>%
    mutate(layer = "total_deaths", count = total_deaths)
  
  red_bar <- plot_data %>%
    mutate(layer = paste0("res_", pathogen_name), count = resistant_deaths)
  
  plot_df <- bind_rows(
    gray_bar %>% select(infection_type, count, layer),
    orange_bar %>% select(infection_type, count, layer),
    red_bar %>% select(infection_type, count, layer)
  )
  
  p <- ggplot() +
    geom_col(data = filter(plot_df, layer == "total_patients"),
             aes(y = infection_type, x = count),
             fill = "grey90", width = 0.6, color = NA) +
    geom_col(data = filter(plot_df, layer == "total_deaths"),
             aes(y = infection_type, x = count),
             fill = "grey65", width = 0.4, color = NA) +
    geom_col(data = filter(plot_df, grepl("^res_", layer)),
             aes(y = infection_type, x = count, fill = layer),
             width = 0.2, color = NA) +
    
    scale_fill_manual(
      values = setNames("#cbdde6", paste0("res_", pathogen_name)),
      labels = pathogen_name,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 2500),
      breaks = seq(0, 2500, 500)
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, color = "black", family = "Times New Roman"),
      axis.text.y = element_text(size = 12, color = "black", family = "Times New Roman",
                                 hjust = 1,  margin = margin(r = -10)),
      axis.text.x = element_text(size = 12, color = "black", family = "Times New Roman"),
      legend.position = c(0.8, 0.1),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.box.background = element_rect(color = NA),
      legend.text = element_text(size = 12, color = "black", family = "Times New Roman")
    ) +
    labs(
      x = " ",
      y = NULL,
      title = organism_label  
    )
  
  print(p)
}

p1 <- draw_pathogen_plot("CRA", expression(bolditalic("Acinetobacter")~bold("spp.")))
p3 <- draw_pathogen_plot("CRP", expression(bolditalic("Pseudomonas")~bold("spp.")))
p4 <- draw_pathogen_plot("VRE", expression(bolditalic("Enterococcus")~bold("spp.")))
p5 <- draw_pathogen_plot("MRSA", expression(bolditalic("Staphylococcus aureus")))

# -------
# Add x title
p5 <- p5 +
  labs(x = "Number of patients") +
  theme(
    axis.title.x = element_text(
      size = 12,             
      face = "plain",        
      family = "Times New Roman", 
      margin = margin(t = 6) 
    )
  )


# -------
# Combine
library(gridExtra)
library(grid)

blank <- nullGrob()


layout_mat <- rbind(
  c(1, 2),
  c(3, 4),
  c(5, 6)
)

#
combined_plot <- grid.arrange(
  grobs = list(p1, p2, p3, p4, p5, blank),
  layout_matrix = layout_mat,
  widths = unit(c(1, 1),  "null"),
  heights = unit(c(1, 1, 1), "null")
)

# Save
ggsave("output/pdf/combined_resistant_deaths.pdf", plot = combined_plot,
       width = 12, height = 8, device = cairo_pdf)
