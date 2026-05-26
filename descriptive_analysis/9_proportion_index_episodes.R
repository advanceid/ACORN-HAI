# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr, 
                 dplyr,
                 tidyr,
                 ggpubr,
                 scales,
                 extrafont,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data/data_table_index_new.RData")
baseline <- readRDS("data/clean_data/data_table_index_new.RData")
org <- readRDS("data/clean_data/ast_all_index.RData")

#Load fonts
loadfonts()

#
baseline[baseline == ""] <- NA
org[org == ""] <- NA

###
df_baseline <- baseline[, c("recordid", "infection_types", "country")]
df_baseline$country <- as.character(df_baseline$country)
df_baseline$country[df_baseline$country == "Hong Kong SAR China"] <- "Hong Kong"
df_baseline$country <- as.factor(df_baseline$country)

df_org <- org %>% 
  filter(!is.na(org_new)) %>%
  group_by(recordid)  %>%
  distinct(recordid, org_new)%>% 
  dplyr::summarize(org_combined_new = paste(org_new, collapse = ", "), .groups = 'drop') %>%
  as.data.frame()

#
df <- left_join(df_baseline, df_org, by = "recordid") %>% 
  filter(!is.na(infection_types) & !is.na(org_combined_new))

# Counts by contry and infection types
fre <- df %>%
  group_by(country, infection_types) %>%
  dplyr::summarise(counts_episode = n(), .groups = 'drop') %>%
  left_join(
    df %>% 
      group_by(country) %>% 
      dplyr::summarise(counts_episode_sum = n(), .groups = 'drop'),
    by = "country"
  ) %>%
  mutate(
    proportion = counts_episode / counts_episode_sum,
    country_with_episode = paste(country, " (", counts_episode_sum, ")", sep = "")
  )


# Plot
fre$infection_types <- factor(fre$infection_types, 
                              levels = c("Healthcare-associated BSI", 
                                         "Hospital-acquired BSI", "VAP"))

p <- ggplot(fre, aes(x = proportion, 
                     y = reorder(country_with_episode, proportion), 
                     fill = infection_types)) +
  geom_col(position = position_dodge(width = 0.6), 
           color = "gray30", 
           width = 0.6, 
           linewidth = 0.2) +
  scale_x_continuous(labels = percent_format(), expand = c(0, 0),
                     limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(
    values = c("VAP" = "#39ad76",
               "Hospital-acquired BSI" = "#352245",
               "Healthcare-associated BSI" = "#dada2b"),
    breaks = c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = "black", family = "Times New Roman"),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.key.height = unit(0.6, "cm"),       
    legend.spacing.y = unit(0.4, "cm"),        
    axis.text = element_text(size = 14, color = "black", family = "Times New Roman"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              family = "Times New Roman", margin = margin(b = 7)),
    plot.margin = unit(c(0.2, 0.7, 0.2, 0.2), "cm") 
  ) +
  labs(title = "Proportion of infection syndromes (%)")

# Save figure
CairoPDF(file = "output/figure/proportion_index_episode.pdf", 
         width = 6, height = 13)
print(p)
dev.off()

###
# Patients number
patients_number <- fre %>% 
  group_by(infection_types) %>% 
  dplyr::summarise(total_counts = sum(counts_episode, na.rm = TRUE))

# Organism counts
fre_org <- org %>% 
  filter(!is.na(pathogen_group_combined)) %>%
  group_by(recordid)  %>%
  distinct(recordid, pathogen_group_combined)%>%
  as.data.frame()

#
df_fre_org <- left_join(df_baseline, fre_org, by = "recordid") %>% 
  filter(!is.na(infection_types) & !is.na(pathogen_group_combined)) %>%
  separate_rows(pathogen_group_combined, sep = ",\\s*") 

# GNB, GPB, Fungi summary
organism_counts <- df_fre_org %>%
  group_by(infection_types, pathogen_group_combined) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(pathogen_group_combined %in% c("GNB", "GPB", "Fungi")) %>%
  pivot_wider(
    names_from = pathogen_group_combined,
    values_from = count,
    values_fill = 0
  )

# Save
saveRDS(patients_number, file = "data/clean_data/patients_number.RData")
saveRDS(organism_counts, file = "data/clean_data/organism_counts.RData")
###
