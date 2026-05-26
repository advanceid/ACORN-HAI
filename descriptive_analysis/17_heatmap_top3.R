# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr, 
                 dplyr,
                 tidyr,
                 ggplot2,
                 rworldmap,
                 extrafont,
                 ggtext,
                 grid,
                 cowplot,
                 patchwork,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df_site <- readRDS("data/clean_data/data_table_index_new.RData")
df_ast <- readRDS("data/clean_data/ast_all_index.RData")

# Delete not recorded in baseline
df_ast_ori <- df_ast[which(df_ast$recordid %in% df_site$recordid),]

# Delete not identified organisms
df_org <- df_ast_ori[-which(is.na(df_ast_ori$org_names_all)),]

# Combine data
df <- left_join(df_org[, c("recordid", "org_names_all", 
                           "ris_Carbapenem", "ris_Third-generation cephalosporin")], 
                df_site[, c("recordid", "country", "infection_types")], 
                by = "recordid")

# Load fonts
loadfonts()

# All countries names
all_countries_names <- unique(df$country)

# Extract the correct order of levels from the original data frame
correct_order <- levels(df$infection_types)

# Pre data for plot
process_data <- function(df, org_filter, ris_var) {

  df_filtered <- df %>%
    filter(org_names_all %in% org_filter) %>%
    group_by(country, infection_types) %>%
    dplyr::summarise(counts_org = n(), .groups = 'drop') %>%
    left_join(
      df %>%
        filter(org_names_all %in% org_filter & .data[[ris_var]] == 1) %>%
        group_by(country, infection_types) %>%
        dplyr::summarise(counts_org_ris = n(), .groups = 'drop'),
      by = c("country", "infection_types")
    ) %>%
    mutate(counts_org_ris = ifelse(is.na(counts_org_ris), 0, counts_org_ris),
           proportion_ris = counts_org_ris / counts_org) %>%
    rename(region = country) %>% 
    # Use complete to ensure all combinations of country and infection_types
    complete(region = all_countries_names, 
             infection_types = levels(factor(df$infection_types)), 
             fill = list(counts_org = 0, counts_org_ris = 0, proportion_ris = 0))
  
  # Ensure correct order of factor levels
  df_filtered$infection_types <- factor(df_filtered$infection_types, levels = correct_order)
  
  # Replace 0 with 0.1 
  df_filtered <- df_filtered %>%
    mutate(proportion_ris = ifelse(proportion_ris == 0, 0.1, proportion_ris))
  
  # Cut off proportions
  df_filtered$proportion_ris_cut <- cut(df_filtered$proportion_ris,
                                        breaks = seq(0, 1.0, 0.1),
                                        labels = c("0 to <10%", "10 to <20%",
                                                   "20 to <30%", "30 to <40%", 
                                                   "40 to <50%", "50 to <60%", 
                                                   "60 to <70%", "70 to <80%", 
                                                   "80 to <90%", "90 to ≤100%"),
                                        right = TRUE, 
                                        ordered_result = TRUE)
  
  df_filtered %>%
    split(.$infection_types)
}

# Process data and split by infection_types
ent_spec_org <- c(
  "E. coli", "K. pneumoniae", "Klebsiella", 
  "Enterobacter", "Serratia", "Proteus", "Morganella",
  "Citrobacter", "Providencia", "Raoultella", "Kluyvera",
  "Cronobacter", "Leclercia", "Pluralibacter", 
  "Edwardsiella", "Pantoea", "Escherichia",
  "Kalamiella"
)

df_used_aci <- process_data(df, c("Acinetobacter spp."), "ris_Carbapenem")
df_used_ent_thir <- process_data(df, ent_spec_org, "ris_Third-generation cephalosporin")
df_used_ent_car <- process_data(df, ent_spec_org, "ris_Carbapenem")

# Combine all lists into one
df_used <- c(df_used_aci, df_used_ent_thir, df_used_ent_car)

# Pre data
df_map <- map_data("world")

# Revise "Hong Kong" to "Hong Kong SAR China"
df_map$region[df_map$subregion == "Hong Kong"] <- "Hong Kong SAR China"

# Plot
df_data <- centers <- df_countries <- list()

system.time({
  for (i in seq_along(df_used)) {
    
    df_countries[[i]] <- df_map %>% filter(region %in% df_used[[i]]$region)
    
    centers[[i]] <- df_countries[[i]] %>%
      group_by(region) %>%
      dplyr::summarize(long = mean(range(long)), lat = mean(range(lat)))
    
    df_data[[i]] <- left_join(df_used[[i]], centers[[i]], by = "region")
    
  }
})


# Initialize a list to store the plots
p_list <- list()

# Define a common color palette and levels for `proportion_ris_cut`
col <- c("#66FF66", "#CCFFCC", "#FFFF99", "#FFFF00",
         "#FFCC00", "#FF9900", "#FF3300", "#8B0000",
         "#663333", "#330000")
names(col) <- c("0 to <10%", "10 to <20%",
                "20 to <30%", "30 to <40%", 
                "40 to <50%", "50 to <60%", 
                "60 to <70%", "70 to <80%", 
                "80 to <90%", "90 to ≤100%")

# Loop through the infection types
for (i in 1:length(df_used)) {
  df_plot <- df_used[[i]]
  
  # Extract infection type for the title
  infection_type_title <- unique(df_plot$infection_types)
  total_index_episode <- sum(df_plot$counts_org)
  
  # Merge with world map data
  df_plot <- df_plot %>%
    left_join(df_map, by = "region")
  
  # Base map plot
  p <- ggplot() +
    geom_polygon(data = df_map, aes(long, lat, group = group), 
                 fill = "white", color = "grey85") +
    geom_polygon(data = df_plot, aes(long, lat, group = group, 
                                     fill = proportion_ris_cut), color = "grey60") +
    coord_cartesian(xlim = c(58, 148), ylim = c(-10, 53)) +
    scale_fill_manual(values = col, drop = FALSE) +  
    guides(fill = guide_legend(override.aes = list(colour = NA))) +
    theme_void() + 
    theme(
      panel.background = element_rect(fill = "grey95", color = NA),
      panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.title = element_text(size = 32, face = "bold", 
                                hjust = 0.5, family = "Times New Roman"),
      plot.subtitle = element_text(hjust = 0.5, size = 32,
                                   family = "Times New Roman")
    )
  
  # Conditionally add the title only for the first three plots
  if (i <= 3) {
    p <- p +
      labs(title = infection_type_title)
  }
  
  # Add the subtitle to all plots
  p <- p + labs(subtitle = paste("n =", format(total_index_episode, big.mark = "", scientific = FALSE)))
  
  # Check if "Singapore" is in the data and add zoom-in with annotation
  if ("Singapore" %in% df_plot$region) {
    p <- p + 
      geom_rect(aes(xmin = 103.5, xmax = 104, ymin = 1.26, ymax = 1.45), 
                color = "black", fill = NA, linewidth = 1) +  
      geom_segment(aes(x = 103.5, y = 1.35, xend = 70.5, yend = -5), 
                   color = "black", linewidth = 1, 
                   arrow = grid::arrow(type = "closed", 
                                       length = unit(0.3, "cm"), angle = 15)) +  
      geom_rect(aes(xmin = 59.7, xmax = 70.2, ymin = -7.5, ymax = -2.6), 
                color = "black", fill = "white", linewidth = 1) +
      annotate("text", x = 65, y = 0.4, label = "Singapore", 
               size = 10, family = "Times New Roman", hjust = 0.5)
    
    # Inset plot for zoom-in
    inset_plot <- ggplot() +
      geom_polygon(data = df_plot, aes(long, lat, group = group, fill = proportion_ris_cut), color = "grey50") +
      coord_cartesian(xlim = c(103.65, 104), ylim = c(1.26, 1.45)) +
      scale_fill_manual(values = col, drop = FALSE) +  
      theme_void() +
      theme(legend.position = "none") 
    
    # Add inset to the main plot
    p <- p + annotation_custom(
      grob = ggplotGrob(inset_plot), 
      xmin = 60, xmax = 70, ymin = -7, ymax = -3
    )
  }
  
  # Check if "Brunei" is in the data and add zoom-in with annotation
  if ("Brunei" %in% df_plot$region) {
    p <- p + 
      geom_rect(aes(xmin = 114, xmax = 115.4, ymin = 4, ymax = 5.1), 
                color = "black", fill = NA, linewidth = 1) +  
      geom_segment(aes(x = 115.4, y = 4.5, xend = 134.7, yend = 7.5), 
                   color = "black", linewidth = 1, 
                   arrow = grid::arrow(type = "closed", 
                                       length = unit(0.3, "cm"), angle = 15)) +  
      geom_rect(aes(xmin = 134.7, xmax = 145.4, ymin = 5.5, ymax = 10.3), 
                color = "black", fill = NA, linewidth = 1) +
      annotate("text", x = 140, y = 12.4, label = "Brunei", 
               size = 10, family = "Times New Roman", hjust = 0.5)
    
    # Inset plot for zoom-in
    inset_plot <- ggplot() +
      geom_polygon(data = df_plot, aes(long, lat, group = group, fill = proportion_ris_cut), 
                   color = "grey50") +
      coord_cartesian(xlim = c(114.3, 115.4), ylim = c(4, 5.1)) +
      scale_fill_manual(values = col, drop = FALSE) +  
      theme_void() +
      theme(legend.position = "none")
    
    
    # Add inset to the main plot
    p <- p + annotation_custom(
      grob = ggplotGrob(inset_plot), 
      xmin = 135, xmax = 145, ymin = 6, ymax = 10
    )
  }
  
  # Check if "Hong Kong" is in the data and add zoom-in with annotation
  if ("Hong Kong SAR China" %in% df_plot$region) {
    p <- p + 
      geom_rect(aes(xmin = 113.5, xmax = 114.4, ymin = 22.1, ymax = 22.6), 
                color = "black", fill = NA, linewidth = 1) +  
      geom_segment(aes(x = 114.4, y = 22.3, xend = 134.7, yend = 18), 
                   color = "black", linewidth = 1, 
                   arrow = grid::arrow(type = "closed", 
                                       length = unit(0.3, "cm"), angle = 15)) +  
      geom_rect(aes(xmin = 134.7, xmax = 145.4, ymin = 16.5, ymax = 19.9), 
                color = "black", fill = NA, linewidth = 1) +
      annotate("text", x = 140, y = 22.7, label = "Hong Kong", 
               size = 10, family = "Times New Roman", hjust = 0.5)
    
    # Inset plot for zoom-in
    inset_plot <- ggplot() +
      geom_polygon(data = df_plot, aes(long, lat, group = group, fill = proportion_ris_cut), 
                   color = "grey50") +
      coord_cartesian(xlim = c(113.8, 114.4), ylim = c(22.2, 22.6)) +
      scale_fill_manual(values = col, drop = FALSE) +  
      theme_void() +
      theme(legend.position = "none")
    
    
    # Add inset to the main plot
    p <- p + annotation_custom(
      grob = ggplotGrob(inset_plot), 
      xmin = 135, xmax = 145, ymin = 16.5, ymax = 19.8
    )
  }
  
  
  # Add the plot to the list
  p_list[[i]] <- p
}

# Create a dummy plot for the legend
dummy_data <- data.frame(
  anti_group = factor(rep("Dummy", 10), levels = "Dummy"),
  organism_names = factor(names(col), levels = names(col)),
  per_cut = factor(names(col), levels = names(col)),
  sum_count = rep("", 10)
)

# Create a dummy plot for the legend
legend_plot <- ggplot(dummy_data, aes(anti_group, organism_names)) +
  geom_tile(aes(fill = per_cut)) +
  scale_fill_manual(values = col, drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 32, hjust = 0, vjust = 0.5,
                                margin = margin(t = 0, r = 25, b = 0, l = 0)),
    legend.text = element_text(size = 32, 
                               margin = margin(t = 0, r = 12, b = 0, l = 3, unit = "pt")),
    legend.margin = margin(t = 0, r = 25, b = 0, l = 120, unit = "pt"),
    legend.background = element_rect(fill = "white", color = NA),
    text = element_text(family = "Times New Roman"),
    legend.spacing.x = unit(40, 'pt'),
    legend.key.height = unit(28, "pt"), 
    legend.key.width = unit(40, "pt")  
  ) +
  guides(fill = guide_legend(title = "Proportion of resistant isolates", 
                             nrow = 2))

# Function to extract the legend
get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract the custom legend
legend <- get_legend(legend_plot)

# Create labels as grobs with Times New Roman font
label_list <- list(
  grid.text("CRA", gp = gpar(fontsize = 32, fontface = "bold", fontfamily = "Times New Roman")),
  grid.text("3GCRE", gp = gpar(fontsize = 32, fontface = "bold", fontfamily = "Times New Roman")),
  grid.text("CRE", gp = gpar(fontsize = 32, fontface = "bold", fontfamily = "Times New Roman"))
)

# Combine the labels with the plot list
combined_plot <- wrap_plots(
  wrap_plots(label_list[[1]], p_list[[1]], p_list[[2]], p_list[[3]], ncol = 4, widths = c(0.2, 1, 1, 1)),
  wrap_plots(label_list[[2]], p_list[[4]], p_list[[5]], p_list[[6]], ncol = 4, widths = c(0.2, 1, 1, 1)),
  wrap_plots(label_list[[3]], p_list[[7]], p_list[[8]], p_list[[9]], ncol = 4, widths = c(0.2, 1, 1, 1)),
  ncol = 1
) +
  plot_layout(guides = 'collect') 

# Use patchwork to combine the plots with the legend
final_plot <- combined_plot / legend + 
  plot_layout(heights = c(1, 1, 1, 0.15)) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(5, 5, 0, 0, unit = "pt"))

print(final_plot)

# Save
CairoPDF(file = "output/figure/proportion_resist_top3_h.pdf", 
         width = 27, height = 21)
plot(final_plot)
dev.off()
###
