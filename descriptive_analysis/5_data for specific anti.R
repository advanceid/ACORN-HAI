# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(readxl,
                 magrittr, 
                 dplyr,
                 tidyr,
                 stringr)
})

## Define working directory
wd <- "./"
setwd(wd)

## Load data
baseline <- readRDS("data/clean_data/data_table_index_new.RData")
ast <- readRDS("data/clean_data/ast_all_index.RData")

# Delete not recorded in baseline
df_ast_ori <- ast[which(ast$recordid %in% baseline$recordid),]

# Delete not identified organisms
ast_used <- df_ast_ori[-which(is.na(df_ast_ori$org_names_all)),]

# Anti group (UNIFIED: expect columns 'anti_group' and 'anti_names')
anti <- read_excel("data/anti group/anti_group_used.xlsx")
org_anti_class <- read_excel("data/anti group/organism_group.xlsx") %>% 
  as.data.frame()

#####
# Create a simplified mapping of ris_ variable to anti_class
ris_col <- grep("^ris_", names(ast_used), value = TRUE)

# Create result list for each ris_ variable
df_del <- lapply(ris_col, function(ris_var) {
  
  anti_class <- gsub("^ris_", "", ris_var)
  
  antibiotics_to_keep <- anti %>%
    dplyr::filter(anti_group == anti_class) %>%
    dplyr::pull(anti_names) %>%
    unique()
  
  # Keep only antibiotic columns that actually exist in ast_used
  antibiotics_to_keep <- intersect(antibiotics_to_keep, names(ast_used))
  
  # Filter the `ast_used` dataframe to include only the required antibiotics and other columns
  selected_columns <- c(ris_var, antibiotics_to_keep, 
                        "infection_types", "org_new", "recordid", "org", "ast_date", "spec_date")
  filtered_data <- ast_used %>%
    dplyr::select(any_of(selected_columns))
  
  return(filtered_data)
})

# Assign names to the list based on `ris_` variable names
names(df_del) <- ris_col

#
count_list <- list()
system.time({
  for (i in seq_along(df_del)) {
    count_list[[i]] <- df_del[[i]] %>%
      dplyr::group_by(!!sym(ris_col[i]), org_new, infection_types) %>%
      dplyr::summarise(
        across(
          .cols = all_of(names(.)[names(.) %in% anti$anti_names]),
          ~ sum(. == 1, na.rm = TRUE),
          .names = "{.col}_count"
        ), 
        total = sum(c_across(all_of(names(.)[names(.) %in% anti$anti_names]))),
        across(
          all_of(names(.)[names(.) %in% anti$anti_names]),
          ~ sum(.)/total,
          .names = "{.col}_percentage"
        ),
        .groups = "drop"
      )
  }
})

#
count_list_new <- list()
pie_result <- data.frame()

system.time({
  for (i in seq_along(count_list)) {
    
    anti_group <- str_remove(names(count_list[[i]])[1], "ris_")
    
    count_list_new[[i]] <- count_list[[i]] %>% 
      dplyr::mutate(anti_group = anti_group) %>%
      dplyr::select("anti_group", starts_with("ris_"), 
                    "org_new", "infection_types", ends_with("_percentage")) %>%
      dplyr::rename_all(~gsub("_percentage", "", .))
    
    colnames(count_list_new[[i]])[2] <- "ris"
    
    pie_result <- dplyr::bind_rows(pie_result, count_list_new[[i]])
  }
})
#
str(pie_result)
for (i in c(2:4)) {pie_result[,i] = as.factor(pie_result[,i])}

# Set labels
pie_result <- within(pie_result, {c(
  
  ris <- factor(ris,
                levels = c(5:1),
                labels = c("NA", "Unknown (U)", "Susceptible (S)",
                           "Intermediate (I)", "Resistant (R)"))
  
)})
#
#####
# To long
count_list_long <- df_list_long <- list()

system.time({
  for (i in seq_along(count_list)) {
    df_used <- count_list[[i]]
    
    anti_group <- str_remove(names(df_used)[1], "ris_")
    
    count_list_long[[i]] <- df_used %>%
      tidyr::pivot_longer(cols = matches("_percentage"),
                          names_to = "antibiotic",
                          values_to = "percentage") %>%
      dplyr::mutate(antibiotic = str_remove(antibiotic, "_percentage"),
                    anti_group = anti_group) 
    
    df_list_long[[i]] <- count_list_long[[i]] %>% 
      dplyr::select("anti_group", starts_with("ris_"), "org_new",
                    "infection_types", "antibiotic", "percentage")
    
    colnames(df_list_long[[i]])[2] <- "ris"
  }
})
#
list_long <- do.call(rbind, df_list_long) %>% 
  as.data.frame()
for (i in 1:5) {list_long[,i] = as.factor(list_long[,i])}

# Set labels
list_long <- within(list_long, {c(
  
  ris <- factor(ris,
                levels = c(5:1),
                labels = c("NA", "Unknown (U)", "Susceptible (S)",
                           "Intermediate (I)", "Resistant (R)"))
  
)})

# Add organism_group
colnames(org_anti_class)[2] <- "org_new"
pie_result <- dplyr::inner_join(org_anti_class,
                                pie_result, 
                                by = c("org_new", "anti_group"))

list_long <- dplyr::inner_join(list_long,
                               org_anti_class, 
                               by = c("org_new", "anti_group"))

# Save
saveRDS(pie_result, "data/clean_data/pie_specific_anti_ast.RData")
saveRDS(list_long, "data/clean_data/bar_specific_anti_ast.RData")
###

