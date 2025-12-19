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
baseline <- readRDS("data/clean data/data_table_index_new.RData")
ast <- readRDS("data/clean data/ast_all_index.RData")

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

# =============================================
# candida
candida_df <- df_ast[which(df_ast$org_new == "Candida spp."),]

length(unique(candida_df$recordid))

length(unique(df_ast$recordid))
length(which(df_ast$Class == "GNB"))
length(which(df_ast$Class == "GPB"))
length(which(df_ast$Class == "Fungi"))

round(prop.table(table(df_ast$Class)) * 100, 1)

length(which(df_ast$org_new == "Candida spp."))/11217

candida_df_BSI <- candida_df[-which(candida_df$infection_types == "VAP"),]

length(which(candida_df_BSI$ris_Azole == 1))
nrow(candida_df_BSI)
length(which(candida_df_BSI$ris_Echinocandin == 1))
length(which(candida_df_BSI$ris_Azole == 1))/nrow(candida_df_BSI)
length(which(candida_df_BSI$ris_Echinocandin == 1))/nrow(candida_df_BSI)

candida_df_BSI %>%
  group_by(infection_types) %>%
  summarise(
    n = n(),
    azoles_count = sum(ris_Azole == 1, na.rm = TRUE),
    azoles_pct   = round(mean(ris_Azole == 1, na.rm = TRUE) * 100, 1),
    echino_count = sum(ris_Echinocandin == 1, na.rm = TRUE),
    echino_pct   = round(mean(ris_Echinocandin == 1, na.rm = TRUE) * 100, 1)
  ) 

# =============================================
# MDR proportion
# Infection types
infection_levels <- c("VAP", "Hospital-acquired BSI", "Healthcare-associated BSI")

# Compute MDR status using df_ast and group by org_new
df_all <- df_ast %>%
  select(recordid, infection_types, starts_with("ris_"), org_new) %>%
  filter(!is.na(org_new), org_new != "Others") %>%
  mutate(across(starts_with("ris_"), ~ ifelse(. == 1, 1, 0))) %>%
  mutate(
    n_resistant = rowSums(select(., starts_with("ris_"))),
    MDR = n_resistant >= 3
  ) %>%
  distinct(recordid, infection_types, org_new, MDR)

# Loop through each organism (org_new)
organism_list <- unique(df_all$org_new)
mdr_summary_list <- list()

for (org_name in organism_list) {
  
  df_org <- df_all %>%
    filter(org_new == org_name)
  
  # By infection type
  mdr_by_type <- df_org %>%
    group_by(infection_types) %>%
    summarise(
      total = n(),
      mdr_count = sum(MDR),
      .groups = "drop"
    ) %>%
    complete(infection_types = infection_levels, fill = list(total = 0, mdr_count = 0)) %>%
    mutate(
      mdr_percent = ifelse(total > 0, round(100 * mdr_count / total, 1), NA),
      mdr_str = ifelse(total > 0, paste0(mdr_percent, "% (", mdr_count, "/", total, ")"), "NA")
    )
  
  # Overall
  total_all <- nrow(df_org)
  mdr_all <- sum(df_org$MDR)
  mdr_overall_str <- ifelse(
    total_all > 0,
    paste0(round(100 * mdr_all / total_all, 1), "% (", mdr_all, "/", total_all, ")"),
    "NA"
  )
  
  # Output row
  row_out <- tibble(
    organism = org_name,
    MDR_VAP = mdr_by_type$mdr_str[mdr_by_type$infection_types == "VAP"],
    MDR_Hospital_acquired_BSI = mdr_by_type$mdr_str[mdr_by_type$infection_types == "Hospital-acquired BSI"],
    MDR_Healthcare_associated_BSI = mdr_by_type$mdr_str[mdr_by_type$infection_types == "Healthcare-associated BSI"],
    MDR_Overall = mdr_overall_str
  )
  
  mdr_summary_list[[org_name]] <- row_out
}

# Combine output
mdr_summary_df <- bind_rows(mdr_summary_list)

# View result
print(mdr_summary_df)

# =============================================
# The commonest antibiotic class resistant
# Get ris_ variables and antibiotic class names
ris_vars <- names(df)[grepl("^ris_", names(df))]
anti_names <- gsub("^ris_", "", ris_vars)
name_map <- setNames(anti_names, ris_vars)

# Identify MDR cases (resistant to ≥3 classes)
df_mdr <- df %>%
  mutate(
    n_resistant = rowSums(across(all_of(ris_vars))),
    MDR = n_resistant >= 3
  ) %>%
  filter(MDR == TRUE)

# Loop for both stratified and overall summaries
mdr_single_class_list <- list()

# Stratified by organism + infection_type
for (org in unique(df_mdr$org_new)) {
  for (inf_type in unique(df_mdr$infection_types)) {
    
    df_sub <- df_mdr %>%
      filter(org_new == org, infection_types == inf_type)
    
    if (nrow(df_sub) == 0) next
    
    resistance_counts <- colSums(df_sub[ris_vars], na.rm = TRUE)
    max_index <- which.max(resistance_counts)
    
    mdr_single_class_list[[paste(org, inf_type, sep = "_")]] <- tibble(
      organism = org,
      infection_type = inf_type,
      most_common_resistant_class = anti_names[max_index],
      count = resistance_counts[max_index],
      total_MDR_cases = nrow(df_sub),
      proportion_percent = round(100 * resistance_counts[max_index] / nrow(df_sub), 1)
    )
  }
}

# Overall (ignore infection type)
for (org in unique(df_mdr$org_new)) {
  
  df_sub <- df_mdr %>%
    filter(org_new == org)
  
  if (nrow(df_sub) == 0) next
  
  resistance_counts <- colSums(df_sub[ris_vars], na.rm = TRUE)
  max_index <- which.max(resistance_counts)
  
  mdr_single_class_list[[paste0(org, "_overall")]] <- tibble(
    organism = org,
    infection_type = "Overall",
    most_common_resistant_class = anti_names[max_index],
    count = resistance_counts[max_index],
    total_MDR_cases = nrow(df_sub),
    proportion_percent = round(100 * resistance_counts[max_index] / nrow(df_sub), 1)
  )
}

# Combine results
mdr_commonest_class_all <- bind_rows(mdr_single_class_list)

mdr_commonest_class_all_overall <- mdr_commonest_class_all[which(mdr_commonest_class_all$infection_type == "Overall"),] 


# =============================================
# Most common resistant combinations among MDR cases
# Compute number of resistant classes
df_mdr <- df %>%
  mutate(
    n_resistant = rowSums(across(all_of(ris_vars))),
    MDR = n_resistant >= 3
  ) %>%
  filter(MDR == TRUE)

# Loop through organisms
mdr_combination_list <- list()

for (org in unique(df_mdr$org_new)) {
  
  df_org <- df_mdr %>%
    filter(org_new == org)
  
  if (nrow(df_org) == 0) next
  
  combo_strings <- apply(df_org[ris_vars], 1, function(x) {
    resistant_classes <- names(x)[x == 1]
    sorted_class <- sort(name_map[resistant_classes])
    paste(sorted_class, collapse = ", ")
  })
  
  combo_df <- tibble(
    organism = org,
    combination = combo_strings
  ) %>%
    count(combination, sort = TRUE)
  
  top_combo <- combo_df[1, ]
  total_mdr <- sum(combo_df$n)
  
  mdr_combination_list[[org]] <- tibble(
    organism = org,
    most_common_resistant_combination = top_combo$combination,
    count = top_combo$n,
    total_MDR_cases = total_mdr,
    proportion_percent = round(100 * top_combo$n / total_mdr, 1)
  )
}

mdr_commonest_combo_df <- bind_rows(mdr_combination_list)
print(mdr_commonest_combo_df)
