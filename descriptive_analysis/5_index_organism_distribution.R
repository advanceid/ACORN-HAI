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
                 RColorBrewer,
                 ggforce,
                 grid,
                 cowplot,
                 patchwork,
                 extrafont,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
baseline <- readRDS("data/clean data/data_table_index_new.RData")
org <- readRDS("data/clean data/ast_all_index.RData")

#Load fonts
loadfonts()

#
baseline[baseline == ""] <- NA
org[org == ""] <- NA

###
df_baseline <- baseline[, c("recordid", "infection_types", "country")]

df_org <- org %>% 
  filter(!is.na(org_new)) %>%
  group_by(recordid)  %>%
  distinct(recordid, org_new)%>% 
  dplyr::summarize(org_combined_new = paste(org_new, collapse = ", "), .groups = 'drop') %>%
  as.data.frame()

#
df <- left_join(df_baseline, df_org, by = "recordid") %>% 
  filter(!is.na(infection_types) & !is.na(org_combined_new))

# Set factor
df_factor <- c("infection_types", "country", "org_combined_new")

df <- df %>%
  mutate(across(all_of(df_factor), as.factor))
  
df_value <- df %>%
  separate_rows(org_combined_new, sep = ", ") %>%
  group_by(infection_types, country, org_combined_new) %>%
  dplyr::summarise(value = n(), .groups = 'drop') %>% 
  filter(!is.na(infection_types) & !is.na(org_combined_new))

df_sum <- df_value %>%
  group_by(infection_types, country) %>%
  dplyr::summarise(value_sum = sum(value), .groups = 'drop')

result <- left_join(df_value, df_sum, by = c("infection_types", "country")) %>% 
  dplyr::rename(region = country) %>% 
  mutate(value_ratio = value / value_sum) %>% 
  as.data.frame()

# reorder
result$org_combined_new <- factor(result$org_combined_new, 
                                  levels = c("E. coli", "K. pneumoniae", 
                                             "Acinetobacter spp.", "Pseudomonas spp.",
                                             "Serratia/Proteus/Morganella/Enterobacter", 
                                             "Stenotrophomonas spp.", "S. aureus", 
                                             "Enterococcus spp.", "Candida spp.", "Others"))

result$infection_types <- factor(result$infection_types, 
                                 levels = c("VAP", 
                                            "Hospital-acquired BSI", 
                                            "Healthcare-associated BSI"))

result$org_combined_fmt <- case_when(
  # Genus + "spp." or "sp." → italicize genus, keep spp. in plain
  grepl("spp\\.$|sp\\.$", result$org_combined_new) ~ 
    gsub("^(\\S+)\\s(spp\\.|sp\\.)$", "italic('\\1')~plain(' \\2')", result$org_combined_new),
  
  # Genus + species (e.g., E. coli) → italicize both, keep space
  grepl("^\\S+\\s\\S+$", result$org_combined_new) ~ 
    gsub("^(\\S+)\\s(\\S+)$", "italic('\\1 \\2')", result$org_combined_new),
  
  # Keep "Others" as plain text
  result$org_combined_new %in% c("Others") ~ result$org_combined_new,
  
  # Default: italicize single-word genus names
  TRUE ~ paste0("italic('", result$org_combined_new, "')")
)

# Pre
df_map <- map_data("world")

# Revise "Hong Kong" to "Hong Kong SAR China"
df_map$region[df_map$subregion == "Hong Kong"] <- "Hong Kong SAR China"

###
df_used <- split(result, result$infection_types)

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

# Set colors
colors <- c("E. coli" = "#ea5c6f",
            "K. pneumoniae" = "#f7905a",
            "Acinetobacter spp." = "#e187cb",
            "Pseudomonas spp." = "#e2b159",
            "Serratia/Proteus/Morganella/Enterobacter" = "#ebed6f",
            "Stenotrophomonas spp." = "#7ee7bb",
            "S. aureus" = "#cee6c8",
            "Enterococcus spp." = "#a9dce6",
            "Candida spp." = "#e4b7d6",
            "Others" = "#9a9a9a")


# Function to extract the legend safely
get_legend <- function(myggplot, return_all = FALSE) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  guides <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (return_all) {
    return(tmp$grobs[guides])
  } else {
    return(tmp$grobs[[guides[1]]])
  }
}


# Construct color mapping with formatted names, preserving original order
label_lookup <- result %>%
  select(org_combined_new, org_combined_fmt) %>%
  distinct() %>%
  filter(org_combined_new %in% names(colors))

# Get formatted names in original color order
fmt_names <- label_lookup$org_combined_fmt[match(names(colors), label_lookup$org_combined_new)]

# Named color vector with formatted labels
colors_fmt <- setNames(colors, fmt_names)

# Ensure the legend respects this order
result$org_combined_fmt <- factor(result$org_combined_fmt, levels = names(colors_fmt))

# Plot
system.time({
  p_list <- list()
  
  for (i in 1:3) {
    
    # Main map without internal title
    map_plot <- ggplot() +
      geom_polygon(data = df_map, aes(long, lat, group = group),
                   fill = "white", color = "grey90") +
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
          fill = org_combined_new
        ),
        color = NA
      ) +
      scale_fill_manual(values = colors, name = " ") +
      theme(legend.position = "none")
    
    # Empty vertical title (for spacing)
    title_grob <- textGrob(
      label = "",  # No label
      rot = 90,
      gp = gpar(fontsize = 14, fontfamily = "Times New Roman"),
      just = "center"
    )
    
    # Combine title and map horizontally
    p_list[[i]] <- plot_grid(title_grob, map_plot,
                             ncol = 2, rel_widths = c(0.04, 1))
  }
  
  # Legend plot to extract
  df_pie[[3]]$org_combined_fmt <- factor(df_pie[[3]]$org_combined_fmt, levels = names(colors_fmt))
  
  legend_plot <- ggplot() + 
    geom_arc_bar(
      data = df_pie[[3]],
      aes(
        x0 = long, y0 = lat,
        r0 = 0, r = 1.75,
        start = 2 * pi * ymin,
        end = 2 * pi * ymax,
        fill = org_combined_fmt
      ),
      color = NA
    ) +
    scale_fill_manual(values = colors_fmt,
                      labels = function(x) parse(text = x),
                      name = "") +
    guides(fill = guide_legend(nrow = 5)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 14, 
                                 margin = margin(r = 1.5, l = 0, unit = "pt")),
      legend.key.size = unit(0.5, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      text = element_text(family = "Times New Roman")
    )
  
  # Get clean legend
  legend <- get_legend(legend_plot)
  
  # Combine maps and legend
  combined_plot_v <- plot_grid(plotlist = p_list, ncol = 1)
  final_plot <- plot_grid(combined_plot_v, legend,
                          ncol = 1, rel_heights = c(10, 1.3))
  
  # Top title: "Organisms"
  main_title <- textGrob(
    label = "Organisms",
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Times New Roman"),
    just = "center"
  )
  
  # Final composition: title + maps + legend
  final_plot_with_title <- plot_grid(
    main_title,
    final_plot,
    ncol = 1,
    rel_heights = c(0.02, 1)
  )
  
  # Draw
  print(final_plot_with_title)
})

# Save
CairoPDF(file = "output/figure/index_organism_distribution_v.pdf", 
         width = 6, height = 13.5)
plot(final_plot_with_title)
dev.off()
###
