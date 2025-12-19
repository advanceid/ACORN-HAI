# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(openxlsx,
                 magrittr, 
                 dplyr)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
baseline <- readRDS("data/clean data/data_table_index_new.RData")
ast <- readRDS("data/clean data/ast_all_index.RData")
org_anti_class <- read.xlsx("data/anti group/organism_group.xlsx") %>% 
  as.data.frame()

# Delete not recorded in baseline
df_ast_ori <- ast[which(ast$recordid %in% baseline$recordid),]

# Delete not identified organisms
micro <- df_ast_ori[-which(is.na(df_ast_ori$org_names_all)),]


length(unique(micro$recordid))

length(unique(df_ast_ori$recordid))

###
options(timeout = 10000)

# Index data
patterns <- c("infection_types", "org_new", "ris_")

df_micro <- micro[, which(grepl(paste(patterns, collapse='|'), names(micro)))] %>%
  as.data.frame()

# 
df_micro <- df_micro[, c(29, 1:28, 30)]

# Remove "ris_" prefix from anti group names
ris_list <- colnames(df_micro[,c(2:29)])
ris_list_clean <- sub("^ris_", "", ris_list)

# Set factor
for (i in c(1:30)) {df_micro[,i] = as.factor(df_micro[,i])}
for (i in c(1:3)) {org_anti_class[,i] = as.factor(org_anti_class[,i])}
###
# Counts based on org_new, antibiotics and infection types
count_list <- list()

system.time({
  for (i in c(2:29)) {
    count_list[[i-1]] <- df_micro %>%
      group_by(df_micro[[1]], df_micro[[i]], df_micro[[30]]) %>%
      dplyr::summarise(count = n())
  }
})

# Create pre data for plot
ad_count_list <- list()

system.time({
  for (i in seq_along(count_list)) {
    count_df <- count_list[[i]]
    count_df$anti_group <- as.factor(rep(i, nrow(count_list[[i]])))
    ad_count_list[[i]] <- count_df
  }
})

#
merged_count <- Reduce(function(x, y) merge(x, y, by = names(ad_count_list[[1]]), all = TRUE), ad_count_list)
#
for (i in c(1:3, 5)) {df_micro[,i] = as.factor(df_micro[,i])}
#
colnames(merged_count)[1:3] <- c("organism_names", 
                                 "anti_susceptibility", 
                                 "infection_types")
# Delete NA in AST result
merged_count <- merged_count[-which(merged_count$anti_susceptibility == 5),]
###
# Pre for proportion of organism
df_types <- split(merged_count, merged_count$infection_types)
df_sum <- all_sum <- df_types_new <- list()

system.time({
  for (i in 1:3) {
    df_sum[[i]] <- df_types[[i]] %>% 
      group_by(organism_names, anti_group) %>%
      dplyr::summarise(df_all_sum = sum(count))
  }
})
#
system.time({
  for (i in 1:3) {
    
    all_sum[[i]] <- matrix(NA, nrow = nrow(df_types[[i]]), ncol = 2)
    colnames(all_sum[[i]]) <- c("sum_count", "per_count")
    
    for (m in 1:nrow(df_sum[[i]])) {
      for (n in 1:nrow(df_types[[i]])) {
        if ((df_types[[i]]$organism_names[n] == df_sum[[i]]$organism_names[m]) & (df_types[[i]]$anti_group[n] == df_sum[[i]]$anti_group[m])) {
          all_sum[[i]][n,1] <- df_sum[[i]]$df_all_sum[m]
          all_sum[[i]][n,2] <- df_types[[i]]$count[n]/all_sum[[i]][n,1]
        } else {
          
        }
      }
    }
  }
})
#
system.time({
  for (i in 1:3) {
    df_types_new[[i]] <- cbind(df_types[[i]], all_sum[[i]])
  }
})
# Set labels
system.time({
  for (i in 1:3) {
    df_types_new[[i]] <- within(df_types_new[[i]], {c(
      
      anti_susceptibility <- factor(anti_susceptibility,
                                    levels = c(4:1),
                                    labels = c("Unknown (U)", 
                                               "Susceptible (S)",
                                               "Intermediate (I)", 
                                               "Resistant (R)")),
      infection_types <- factor(infection_types, 
                                levels = c("VAP",
                                           "Hospital-acquired BSI",
                                           "Healthcare-associated BSI")),
      anti_group <- factor(anti_group, 
                           levels = c(1:28), 
                           labels = ris_list_clean)
      
    )})
  }
})
###
# Delete irrelevant antibiotics
df_types_new_del <- list()

system.time({
  for (i in seq_along(df_types_new)) {
    df_types_new_del[[i]] <- inner_join(df_types_new[[i]],
                                        org_anti_class, 
                                        by = c("organism_names", "anti_group"))
    
  }
})
# Save 
saveRDS(df_types_new_del, file = "data/clean data/percent_RIS.RData")
###