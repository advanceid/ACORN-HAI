# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr, 
                 dplyr,
                 tibble,
                 tidyr,
                 extrafont,
                 ggforce,
                 ggupset,
                 ggplot2,
                 gridExtra,
                 grid,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
baseline <- readRDS("data/clean_data/data_table_index_new.RData")
ast <- readRDS("data/clean_data/ast_all_index.RData")

# Delete not recorded in baseline
df_ast <- ast[which(ast$recordid %in% baseline$recordid),]

# Delete not identified organisms
df_ast <- df_ast[-which(is.na(df_ast$org_names_all)),]
df_ast <- df_ast[df_ast$Class != "Unknown", ]

# Load fonts
loadfonts()

# Infection types
df_ast$infection_types <- factor(df_ast$infection_types,
                                 levels = c("VAP",
                                            "Hospital-acquired BSI",
                                            "Healthcare-associated BSI"))

# Antibiotic classes
df <- df_ast %>%
  select(recordid, infection_types, starts_with("ris_"), org_new) %>%
  na.omit() %>%
  filter(org_new != "Others") %>%
  mutate(across(starts_with("ris_"), ~ ifelse(. == 1, 1, 0))) %>% 
  filter(rowSums(select(., starts_with("ris_"))) > 0) 


# Split data by infection types and then by org_new
df_used <- split(df, list(df$infection_types, df$org_new), drop = TRUE)

# ------------------------------------------------------------------------------
# Check unique(recordid) and nrow
result <- lapply(df_used, function(sub_df) {
  n_unique_recordid <- length(unique(sub_df$recordid))
  n_rows <- nrow(sub_df)
  list(
    n_unique_recordid = n_unique_recordid,
    n_rows = n_rows,
    same = n_unique_recordid == n_rows
  )
})

# 
for (i in seq_along(result)) {
  cat("Dataframe", i, ":\n")
  cat("Unique recordid count:", result[[i]]$n_unique_recordid, "\n")
  cat("Total row count:", result[[i]]$n_rows, "\n")
  cat("Same:", result[[i]]$same, "\n\n")
}

# ------------------------------------------------------------------------------
# Initialize an empty list to store combined results
combined_ris_list <- list()

# Get unique org_new
org_new_list <- unique(unlist(lapply(df_used, function(x) unique(x$org_new))))

# Loop through each unique org_new
for (org_new in org_new_list) {
  # Filter data frames with the current org_new
  filtered_dfs <- lapply(df_used, function(sub_df) {
    if (unique(sub_df$org_new) == org_new) {
      return(sub_df)
    } else {
      return(NULL)
    }
  })
  filtered_dfs <- filtered_dfs[!sapply(filtered_dfs, is.null)]
  
  # Process each filtered data frame
  tidy_sub_dfs <- lapply(filtered_dfs, function(sub_df) {
    # Select the necessary columns
    sub_df <- sub_df %>%
      select(recordid, starts_with("ris_"), infection_types) %>%
      rename_with(~ gsub("^ris_", "", .), starts_with("ris_"))
    
    # Extract the infection_types column from the original dataframe
    infection_type <- unique(sub_df$infection_types)
    
    # Make recordid values unique in the format A1, A2, ...
    sub_df <- sub_df %>%
      mutate(recordid = make.unique(paste0("A", seq_len(nrow(sub_df))))) %>%
      column_to_rownames(var = "recordid") %>% 
      select(-infection_types)
    
    # Transpose the data frame and convert it back to a data frame
    sub_df <- t(sub_df) %>% 
      as.data.frame() %>%
      rownames_to_column("anti_class")
    
    # Remove columns that contain only zeros
    sub_df <- sub_df[, colSums(sub_df != 0) > 0]
    
    # Transform the data to long format
    tidy_sub_df <- sub_df %>%
      pivot_longer(cols = -c(anti_class), names_to = "ID", values_to = "Member") %>%
      filter(Member != 0) %>%
      select(-Member)
    
    # Add the infection_types column back
    tidy_sub_df$infection_types <- infection_type
    
    return(tidy_sub_df)
  })
  
  # Combine all tidy_sub_dfs for the current org_new
  combined_tidy_sub_df <- do.call(rbind, tidy_sub_dfs)
  
  # Store the combined result in the list
  combined_ris_list[[org_new]] <- combined_tidy_sub_df
}

###
# Initialize a list to store data and plots
processed_data_ori_save <- p_list <- list()

# Iterate over combined_ris_list and create plots
for (org_new in names(combined_ris_list)) {
  
  processed_data_ori <- combined_ris_list[[org_new]] %>%
    group_by(ID, infection_types) %>%
    summarise(anti_class = list(anti_class), .groups = 'drop') %>%
    ungroup() %>%
    mutate(anti_class_str = sapply(anti_class, toString)) %>%
    group_by(infection_types, anti_class_str) %>%
    mutate(count = n()) %>%
    arrange(infection_types)
  
  anti_class_list <- unlist(strsplit(processed_data_ori$anti_class_str, ", "))
  anti_class_frequency <- as.data.frame(table(anti_class_list))
  colnames(anti_class_frequency) <- c("anti", "frequency")
  
  frequency_map <- setNames(anti_class_frequency$frequency, anti_class_frequency$anti)
  
  update_anti_class <- function(anti_class_str, frequency_map) {
    classes <- unlist(strsplit(anti_class_str, ", "))
    updated_classes <- sapply(classes, function(class) {
      freq <- frequency_map[[class]]
      paste(class, "(", freq, ")", sep = "")
    })
    paste(updated_classes, collapse = ", ")
  }
  
  processed_data_ori$anti_class_str <- sapply(processed_data_ori$anti_class_str, update_anti_class, frequency_map)
  
  processed_data_ori$anti_class_list <- strsplit(processed_data_ori$anti_class_str, ", ")
  
  total_index_episode <- length(processed_data_ori$ID)
  
  # If org_new is "Stenotrophomonas spp.", remove rows with Carbapenem or Third-generation cephalosporin
  if (org_new == "Stenotrophomonas spp.") {
    processed_data_ori <- processed_data_ori %>%
      filter(!sapply(anti_class_list, function(x) any(grepl("Carbapenem|Third-generation cephalosporin", x))))
  }
  
  
  processed_data_ori_save[[org_new]] <- processed_data_ori
  
  processed_data <- processed_data_ori %>% filter(count >= 3)
  
  max_count <- max(processed_data$count)
  break_interval <- ceiling(max_count / 10)
  
  # Create plot title expression
  plot_title_expr <- if (grepl("spp\\.$", org_new)) {
    genus <- gsub(" spp\\.$", "", org_new)
    bquote(bolditalic(.(genus)) ~ bold(" spp."))
  } else if (org_new == "Others") {
    bquote(bold(.(org_new)))
  } else {
    bquote(bolditalic(.(org_new)))
  }
  
  p <- ggplot(processed_data, aes(x = anti_class_list)) +
    geom_bar(fill = "#0073C2FF", color = NA) +
    scale_x_upset() +
    scale_y_continuous(breaks = seq(0, max_count, by = break_interval),
                       expand = c(0, 0),
                       position = "right") +  # Move y-axis to the right
    facet_grid(infection_types ~ ., scales = "free", space = "free", 
               switch = "y") +  # Move facet labels to the left
    theme_minimal(base_size = 15) +
    theme_combmatrix(combmatrix.label.text = element_text(color = "black", 
                                                          size = 16, family = "Times New Roman"),
                     combmatrix.label.make_space = FALSE,
                     plot.margin = unit(c(2, 3, 2, 188), "pt")) +
    theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
          panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5, face = "bold.italic", 
                                    size = 16, family = "Times New Roman"),
          plot.subtitle = element_text(hjust = 0.5, size = 16, family = "Times New Roman"),
          axis.title.y = element_text(face = "bold", size = 16, family = "Times New Roman"),
          axis.text.y = element_text(color = "black", size = 16, family = "Times New Roman"),
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 16), 
          strip.placement = "outside",
          strip.switch.pad.grid = unit(0.6, "lines"), 
          panel.spacing = unit(1, "lines"),
          text = element_text(family = "Times New Roman")) +
    labs(
      x = "",
      y = "Index episodes",
      title = plot_title_expr,
      subtitle = paste("n =", format(total_index_episode, big.mark = "", scientific = FALSE))
    )
  
  p_list[[org_new]] <- p
  print(p_list[[org_new]])
}

# Save
for (i in c(1:9)) {
  assign(paste0("p", i), p_list[[i]])
}

# Create an empty grob for the first column
empty_col <- nullGrob() 

# Enterobacteriaceae
row1 <- arrangeGrob(p1, p5, ncol = 2, widths = c(1.2, 0.8))
row2 <- arrangeGrob(empty_col, p2, empty_col, ncol = 3, widths = c(0.04, 1, 0.06))
plot_1 <- arrangeGrob(row1, row2, ncol = 1, heights = c(1.1, 1.4))

# Non-fermenters
row3 <- arrangeGrob(empty_col, p3, empty_col, ncol = 3, widths = c(0.07, 1, 0.03))
row4 <- arrangeGrob(empty_col, p4, empty_col, p6, empty_col, ncol = 5, 
                    widths = c(0.05, 0.9, 0.07, 0.23, 0.02))
plot_2 <- arrangeGrob(row3, row4, ncol = 1, heights = c(1.2, 1))

# GPB, fungi
row5 <- arrangeGrob(empty_col, p7, empty_col, ncol = 3, widths = c(0.17, 1, 0.17))
row6 <- arrangeGrob(p8, p9, ncol = 2, 
                    widths = c(0.9, 0.6))
plot_3 <- arrangeGrob(row5, row6, ncol = 1, heights = c(1, 1))


# Save the plot
ggsave("gnb1_MDROS.pdf", 
       plot = plot_1,
       width = 24, height = 18, 
       device = cairo_pdf,
       path = "output/figure/",
       limitsize = FALSE)

ggsave("gnb2_MDROS.pdf", 
       plot = plot_2,
       width = 24, height = 18, 
       device = cairo_pdf,
       path = "output/figure/",
       limitsize = FALSE)

ggsave("gpb_fungi_MDROS.pdf", 
       plot = plot_3,
       width = 17, height = 17, 
       device = cairo_pdf,
       path = "output/figure/",
       limitsize = FALSE)

