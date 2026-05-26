# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr,
                 dplyr,
                 tidyr,
                 ggplot2,
                 extrafont,
                 rworldmap,
                 RColorBrewer,
                 geosphere,
                 ggforce,
                 cowplot,
                 grid,
                 gridExtra,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data/data_table_index_new.RData")

#Load fonts
loadfonts() 

# Delete NA
df[df == ""] <- NA
df <- df[, c("recordid", "infection_types", "country", "pathogen_combined_types")]%>% 
  filter(!is.na(infection_types) & !is.na(pathogen_combined_types))

# Set factor
df_factor <- c("infection_types", "country", "pathogen_combined_types")

df <- df %>%
  mutate(across(all_of(df_factor), as.factor))

df_value <- df %>%
  separate_rows(pathogen_combined_types, sep = ", ") %>%
  group_by(infection_types, country, pathogen_combined_types) %>%
  dplyr::summarise(value = n(), .groups = 'drop') %>% 
  filter(!is.na(infection_types) & !is.na(pathogen_combined_types))

df_sum <- df_value %>%
  group_by(infection_types, country) %>%
  dplyr::summarise(value_sum = sum(value), .groups = 'drop')

result <- left_join(df_value, df_sum, by = c("infection_types", "country")) %>% 
  dplyr::rename(region = country) %>% 
  mutate(value_ratio = value / value_sum) %>% 
  as.data.frame()

# reorder
result$pathogen_combined_types <- factor(result$pathogen_combined_types, 
                                         levels = c("Gram-negative bacteria", 
                                                    "Gram-positive bacteria", 
                                                    "Fungi",
                                                    "Polymicrobial"))

result$infection_types <- factor(result$infection_types, 
                                 levels = c("VAP", 
                                            "Hospital-acquired BSI", 
                                            "Healthcare-associated BSI"))
###
df_map <- map_data("world")

# Revise "Hong Kong" to "Hong Kong SAR China"
df_map$region[df_map$subregion == "Hong Kong"] <- "Hong Kong SAR China"

###
df_used <- split(result, result$infection_types)

# Function to calculate proportions
calculate_proportions <- function(data) {
  pathogen_sum <- aggregate(value ~ pathogen_combined_types, data, sum)
  unique_data <- data[!duplicated(data[c("region", "value_sum")]), ]
  total_value_sum <- sum(unique_data$value_sum)
  pathogen_sum$ratio <- pathogen_sum$value / total_value_sum
  return(pathogen_sum)
}

# Apply the function to each subset
proportion_results <- lapply(df_used, calculate_proportions)

# Display the result
print(proportion_results)

#
df_pie <- df_data <- centers <- df_countries <- list()

# Adjust centers
adjust_centers <- function(df, min_dist_km = 200, jitter_step = 0.3, max_iter = 100) {
  adjusted <- df
  
  for (iter in 1:max_iter) {
    changed <- FALSE

    for (i in 1:(nrow(adjusted) - 1)) {
      for (j in (i + 1):nrow(adjusted)) {
        lon1 <- adjusted$long[i]; lat1 <- adjusted$lat[i]
        lon2 <- adjusted$long[j]; lat2 <- adjusted$lat[j]
        
        if (any(is.na(c(lon1, lat1, lon2, lat2)))) next
        
        dist <- geosphere::distGeo(c(lon1, lat1), c(lon2, lat2)) / 1000
        
        if (dist < min_dist_km) {

          dx <- lon2 - lon1
          dy <- lat2 - lat1
          norm <- sqrt(dx^2 + dy^2)
          if (norm == 0) {
            angle <- runif(1, 0, 2 * pi)
            dx <- cos(angle)
            dy <- sin(angle)
          } else {
            dx <- dx / norm
            dy <- dy / norm
          }
          
          adjusted$long[i] <- adjusted$long[i] - dx * jitter_step
          adjusted$lat[i]  <- adjusted$lat[i] - dy * jitter_step
          adjusted$long[j] <- adjusted$long[j] + dx * jitter_step
          adjusted$lat[j]  <- adjusted$lat[j] + dy * jitter_step
          
          changed <- TRUE
        }
      }
    }
    
    if (!changed) break 
  }
  
  return(adjusted)
}


# Add adjust
system.time({
  for (i in seq_along(df_used)) {
    
    df_countries[[i]] <- df_map %>% filter(region %in% df_used[[i]]$region)
    
    centers_raw <- df_countries[[i]] %>%
      group_by(region) %>%
      dplyr::summarize(long = mean(range(long)), lat = mean(range(lat)))
    
    centers[[i]] <- adjust_centers(centers_raw, min_dist_km = 200)
    
    df_data[[i]] <- left_join(df_used[[i]], centers[[i]], by = "region")
    
    df_pie[[i]] <- df_data[[i]] %>%
      group_by(region, long, lat) %>%
      mutate(
        cum_value = cumsum(value_ratio),
        ymax = cum_value,
        ymin = lag(cum_value, default = 0)
      )
  }
})

# 
system.time({
  for (i in seq_along(df_used)) {
    
    df_countries[[i]] <- df_map %>% filter(region %in% df_used[[i]]$region)
    
    centers[[i]] <- df_countries[[i]] %>%
      group_by(region) %>%
      dplyr::summarize(long = mean(range(long)), lat = mean(range(lat)))
    
    df_data[[i]] <- left_join(df_used[[i]], centers[[i]], by = "region")
    
    df_pie[[i]] <- df_data[[i]] %>%
      group_by(region, long, lat) %>%
      mutate(
        cum_value = cumsum(value_ratio),
        ymax = cum_value,
        ymin = lag(cum_value, default = 0)
      )
  }
})

# Function to extract the legend from a ggplot object
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#
colors <- c("Gram-negative bacteria" = "#1B9E77",
            "Gram-positive bacteria" = "#2B83BA",
            "Fungi" = "#D95F02",
            "Polymicrobial" = "#E6AB02")

# Plot
system.time({
  p_list <- list()
  
  for (i in 1:3) {
    title_label <- levels(df_used[[i]]$infection_types)
    
    # Main map (no title)
    map_plot <- ggplot() +
      geom_polygon(data = df_map, aes(long, lat, group = group), fill = "white", color = "grey90") +
      coord_cartesian(xlim = c(58, 148), ylim = c(-10, 53)) +
      theme_void() + 
      theme(
        panel.background = element_rect(fill = "grey90", color = NA),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white", color = NA),
        text = element_text(family = "Times New Roman"),
        plot.margin = margin(r = 4, b = 5, unit = "pt")
      ) +
      geom_arc_bar(
        data = df_pie[[i]],
        aes(
          x0 = long, y0 = lat,
          r0 = 0, r = 1.75,
          start = 2 * pi * ymin,
          end = 2 * pi * ymax,
          fill = pathogen_combined_types
        ),
        color = NA
      ) +
      scale_fill_manual(values = colors, name = " ") +
      theme(legend.position = "none")
    
    # Left vertical title
    title_grob <- textGrob(
      label = title_label[i],
      rot = 90,
      gp = gpar(fontsize = 14, fontfamily = "Times New Roman"),
      just = "center"
    )
    
    # Combine title and map side by side
    p_list[[i]] <- plot_grid(title_grob, map_plot, 
                             ncol = 2, rel_widths = c(0.04, 1))
  }
  
  # Legend
  legend_plot <- ggplot() + 
    geom_arc_bar(
      data = df_pie[[2]],
      aes(
        x0 = long, y0 = lat,
        r0 = 0, r = 1.75,
        start = 2 * pi * ymin,
        end = 2 * pi * ymax,
        fill = pathogen_combined_types
      ),
      color = NA
    ) +
    scale_fill_manual(values = colors, name = "") +
    guides(fill = guide_legend(nrow = 4)) + 
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14, 
                                     margin = margin(r = 10, l = 0, unit = "pt")),
          legend.key.size = unit(0.5, "cm"),
          legend.spacing.x = unit(0.5, "cm"),
          legend.spacing.y = unit(0.1, "cm"),
          text = element_text(family = "Times New Roman"))
  
  legend <- get_legend(legend_plot)
  
  # Combine maps and legend vertically
  combined_plot_v <- plot_grid(plotlist = p_list, ncol = 1)
  final_plot <- plot_grid(combined_plot_v, legend, 
                          ncol = 1, rel_heights = c(10, 1.3))
  
  # Top main title
  main_title <- textGrob(
    label = "Type of pathogens",
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Times New Roman"),
    just = "center"
  )
  
  # Final layout: title + maps + legend
  final_plot_with_title <- plot_grid(
    main_title,
    final_plot,
    ncol = 1,
    rel_heights = c(0.02, 1)
  )
  
  print(final_plot_with_title)
})

# Save
CairoPDF(file = "output/figure/index_pathogen_types_distribution_v.pdf", 
         width = 6, height = 13.5)
plot(final_plot_with_title)
dev.off()
###
