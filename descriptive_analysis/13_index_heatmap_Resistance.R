# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr, dplyr, viridis, ggplot2, ggh4x, grid,
    cowplot, patchwork, extrafont, forcats, Cairo
  )
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# Load data
df <- readRDS("data/clean data/percent_RIS.RData")
index_number <- readRDS("data/clean data/organism_counts.RData")

# Add subgroup: ensure "Resistant (R)" exists for Stenotrophomonas × Tetracyclines
for (i in seq_along(df)) {
  df_i <- df[[i]]
  
  if (is.factor(df_i$anti_susceptibility)) {
    df_i$anti_susceptibility <- forcats::fct_expand(df_i$anti_susceptibility, "Resistant (R)")
  }
  
  df_sub <- df_i %>%
    dplyr::filter(
      organism_names == "Stenotrophomonas spp.",
      anti_group == "Tetracyclines"
    )
  
  has_R <- any(df_sub$anti_susceptibility == "Resistant (R)", na.rm = TRUE)
  
  if (!has_R && nrow(df_sub) > 0) {
    new_line <- df_sub[1, ]
    new_line$anti_susceptibility <- "Resistant (R)"
    new_line$count     <- 0L
    new_line$sum_count <- 0L
    new_line$per_count <- 0
    df_i <- dplyr::bind_rows(df_i, new_line)
  }
  
  df[[i]] <- df_i
}

# Select R and cut proportion
df_used <- list()

system.time({
  for (i in seq_along(df)) {
    
    df[[i]]$anti_group <- droplevels(df[[i]]$anti_group)
    
    levels(df[[i]]$anti_group)[c(2, 7, 10, 13, 15)] <-
      c("Anti-pseudomonal penicillin/\nBeta-lactamase inhibitor",
        "Fourth-generation\ncephalosporin",
        "Other Beta-lactam/\nBeta-lactamase inhibitor",
        "Sulfonamide-\ntrimethoprim-combination",
        "Third-generation\ncephalosporin")
    
    df[[i]] <- df[[i]] %>%
      dplyr::mutate(
        per_count = as.numeric(per_count),
        per_count = dplyr::if_else(is.na(per_count) & sum_count > 0, 0, per_count),
        per_count = pmin(pmax(per_count, 0), 1)
      )
    
    r_row <- which(df[[i]]$anti_susceptibility == "Resistant (R)")
    
    df_used[[i]] <- df[[i]][r_row, ]
    
    df_used[[i]]$per_cut <- cut(
      df_used[[i]]$per_count,
      breaks = c(seq(-0.0001, 1.0, 0.1)),
      levels = 1:10,
      labels = c("0 to ≤10%", "10 to ≤20%",
                 "20 to ≤30%", "30 to ≤40%",
                 "40 to ≤50%", "50 to ≤60%",
                 "60 to ≤70%", "70 to ≤80%",
                 "80 to ≤90%", "90 to ≤100%"),
      right = TRUE,
      ordered_result = TRUE
    )
    
    df_used[[i]]$organism_group <- factor(
      df_used[[i]]$organism_group,
      levels = c("GNB", "GPB", "Fungi")
    )
    
    # Move new order
    new_order <- c("Anti-pseudomonal penicillin/\nBeta-lactamase inhibitor",
                   "Other Beta-lactam/\nBeta-lactamase inhibitor",
                   "Third-generation\ncephalosporin",
                   "Tetracycline", "Fluoroquinolone",
                   "Sulfonamide-\ntrimethoprim-combination",
                   "Glycopeptide", "Macrolide", "Penicillin",
                   "Azole", "Echinocandin", "Polyene")
    
    anti_index <- which(levels(df_used[[i]]$anti_group) %in% new_order)
    levels_reordered <- c(levels(df_used[[i]]$anti_group)[-anti_index], new_order)
    
    # Update levels for anti_group
    df_used[[i]]$anti_group <- factor(df_used[[i]]$anti_group, levels = levels_reordered)
    
    df_used[[i]]$per_cut <- ifelse(is.na(df_used[[i]]$per_cut),
                                   "0 to ≤10%", as.character(df_used[[i]]$per_cut))
  }
})

# Define the color palette and common legend
col <- viridis_pal(option = "plasma", end = 0.9, direction = -1)(10)
names(col) <- c("0 to ≤10%", "10 to ≤20%", "20 to ≤30%", "30 to ≤40%",
                "40 to ≤50%", "50 to ≤60%", "60 to ≤70%", "70 to ≤80%",
                "80 to ≤90%", "90 to ≤100%")

for (i in seq_along(df_used)) {
  if (!is.null(df_used[[i]]) && nrow(df_used[[i]]) > 0) {
    df_used[[i]]$per_cut <- factor(df_used[[i]]$per_cut, levels = names(col), ordered = TRUE)
  }
}

# Create plots
dev.new()

# Function to create plots without individual legends
create_plot <- function(data, title, count, show_x_axis = FALSE, caption = "") {
  
  # Prepare formatted organism labels
  org_labels <- levels(data$organism_names)
  
  org_fmt <- sapply(org_labels, function(x) {
    if (x == "Serratia/Proteus/Morganella/Enterobacter") {
      return(bquote(
        atop(italic("Serratia/Proteus/") * phantom("Q"),
             italic("Morganella/Enterobacter") * vphantom("Q"))
      ))
    } else if (x == "Others") {
      return(x)
    } else if (grepl("spp\\.$", x)) {
      genus <- gsub(" spp\\.$", "", x)
      return(bquote(italic(.(genus)) ~ plain(" spp.")))
    } else if (grepl("/", x)) {
      return(bquote(italic(.(x))))
    } else {
      return(bquote(italic(.(x))))
    }
  })
  names(org_fmt) <- org_labels
  
  sum_data <- data %>%
    dplyr::distinct(organism_names, anti_group, sum_count)
  
  # Format count with thousands separator
  formatted_count <- format(count, big.mark = "", scientific = FALSE)
  plot_title <- paste0(title, " (n = ", formatted_count, ")")
  
  plot <- ggplot(data, aes(anti_group, organism_names)) +
    geom_tile(aes(fill = per_cut)) +
    geom_text(aes(label = sum_count), color = "white", size = 4) +
    scale_fill_manual(values = col, drop = FALSE, na.value = col["0 to ≤10%"]) +
    theme_bw() +
    labs(x = NULL, y = NULL, caption = caption) +
    ggtitle(plot_title) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = if (show_x_axis) element_line(color = "black", linewidth = 0.5) else element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.text.x = if (show_x_axis) element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black", margin = margin(t = 8)) else element_blank(),
      axis.text.y = element_text(color = "black", size = 14),
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      axis.ticks.y = element_line(),
      legend.position = "none",
      strip.text = element_blank(),
      strip.background = element_rect(color = "black", fill = "grey80", linewidth = NA),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 10),
      text = element_text(family = "Times New Roman")
    ) +
    scale_y_discrete(labels = org_fmt) +
    facet_grid(organism_group ~ ., scales = "free_y", space = "free")
  
  return(plot)
}

# Create individual plots
index_number$total_counts <- rowSums(index_number[, 2:4])

p1 <- create_plot(
  df_used[[3]], "VAP",
  index_number$total_counts[index_number$infection_types == "VAP"],
  show_x_axis = FALSE
)
p2 <- create_plot(
  df_used[[2]], "Hospital-acquired BSI",
  index_number$total_counts[index_number$infection_types == "Hospital-acquired BSI"],
  show_x_axis = FALSE
)
p3 <- create_plot(
  df_used[[1]], "Healthcare-associated BSI",
  index_number$total_counts[index_number$infection_types == "Healthcare-associated BSI"],
  show_x_axis = TRUE
)

combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 1)


# Create a dummy plot for the legend
dummy_data <- data.frame(
  anti_group = factor(rep("Dummy", 10), levels = "Dummy"),
  organism_names = factor(names(col), levels = names(col)),
  per_cut = factor(names(col), levels = names(col)),
  sum_count = rep("", 10)
)

legend_plot <- ggplot(dummy_data, aes(anti_group, organism_names)) +
  geom_tile(aes(fill = per_cut)) +
  scale_fill_manual(values = col, drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 14, hjust = 0, vjust = 0.5,
                                margin = margin(t = 0, r = 9, b = 0, l = 2)),
    legend.text = element_text(size = 14, margin = margin(r = 0, l = 1, unit = "pt")),
    legend.margin = margin(b = 5),
    legend.background = element_rect(fill = "white", color = NA),
    text = element_text(family = "Times New Roman"),
    legend.spacing.x = unit(16, 'pt')
  ) +
  guides(fill = guide_legend(title = "Proportion of resistant isolates", 
                             nrow = 2))

# Extract the legend
get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(legend_plot)

# Combine with the legend as another row
final_plot <- combined_plot / legend +
  plot_layout(ncol = 1,  
              heights = c(2.7, 3, 3, 0.3)) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

# Save
CairoPDF("output/figure/Heatmap_R_index.pdf", width = 11, height = 18)
print(final_plot)
dev.off()
###
