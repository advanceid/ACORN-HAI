# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr, 
                 dplyr,
                 ggplot2,
                 forcats,
                 grid,
                 gridExtra,
                 patchwork,
                 extrafont,
                 Cairo)
})
###
# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data/percent_RIS.RData")
index_number <- readRDS("data/clean_data/organism_counts.RData")

# Load fonts
loadfonts()

# Make plots
options(timeout = 60)

data_list <- plot_list <- list()

# Split data
system.time({
  for (i in seq_along(df)) {
    
    df_used <- df[[i]]
    df_used <- df[[i]][which(df[[i]]$organism_group %in% c("GPB", "Fungi")),]
    df_used$organism_names <- droplevels(df_used$organism_names)
    df_used$anti_group <- droplevels(df_used$anti_group)
    
    levels(df_used$anti_group)[8] <- "Sulfonamide-\ntrimethoprim-\ncombination"
    
    # Save
    data_list[[i]] <- df_used
    
    # Move new order
    new_order <- c("Fluoroquinolone", 
                   "Sulfonamide-\ntrimethoprim-\ncombination", "Glycopeptide",
                   "Macrolide", "Penicillin", 
                   "Azole", "Echinocandin", "Polyene")
    
    # Update levels for anti_group and org
    data_list[[i]]$anti_group <- factor(data_list[[i]]$anti_group, levels = new_order)
    
  }
})

temp <- data_list[[1]]
data_list[[1]] <- data_list[[3]]
data_list[[3]] <- temp

index_number$total_counts <- c(index_number$Fungi + index_number$GPB)

# Plot
options(expressions = 50000) 

system.time({
  for (i in 1:3) {
    df_used <- data_list[[i]]
    
    # Generate display labels
    org_labels <- as.character(df_used$organism_names)
    org_fmt <- sapply(org_labels, function(x) {
      if (grepl("spp\\.$", x)) {
        genus <- gsub(" spp\\.$", "", x)
        return(bquote(italic(.(genus)) ~ plain(" spp.")))
      } else if (grepl("/", x)) {
        return(bquote(italic(.(x))))
      } else {
        return(bquote(italic(.(x))))
      }
    })
    names(org_fmt) <- org_labels
    
    tit <- levels(df_used$infection_types)
    total_count <- index_number$total_counts[i]
    
    df_used$organism_names <- factor(df_used$organism_names, 
                                     levels = c("Candida spp.", "S. aureus", "Enterococcus spp."))
    
    # Plot
    p <- ggplot(df_used, 
                aes(y = organism_names,
                    x = per_count, 
                    fill = anti_susceptibility)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      geom_text(aes(label = sum_count, x = 0.9, y = organism_names), 
                vjust = 0.5, size = 4.5, fontface = "bold", color = "white") +
      facet_wrap( ~ anti_group, nrow = 1) +
      scale_fill_manual(
        values = c("Resistant (R)" = "#352245",  
                   "Intermediate (I)" = "#efe13a",
                   "Susceptible (S)" = "#39ad76",
                   "Unknown (U)" = "#9a9a9a"),
        breaks = c("Resistant (R)", "Intermediate (I)", "Susceptible (S)", "Unknown (U)")) +
      scale_x_continuous(labels = scales::percent) +
      theme_bw() + 
      labs(x = "", y = "", title = "") +
      theme(
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 18, hjust = 0, vjust = 0.5,
                                    margin = margin(r = 20, unit = "pt")),
        legend.text = element_text(size = 18, 
                                   margin = margin(l = 3, r = 6, unit = "pt")),
        legend.position = ifelse(i == 3, "bottom", "none"),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = if (i == 3) {
          element_text(angle = 90, 
                       margin = margin(t = 18),
                       hjust = 1, vjust = 0.5, 
                       colour = "black", size = 18)
        } else {
          element_blank()
        },  
        axis.text.y = element_text(margin = margin(r = 2),
                                   hjust = 1, vjust = 0.5, 
                                   colour = "black", size = 18),
        axis.line = element_line(colour = "black", linewidth = 0.45),
        panel.spacing.x = unit(0.5, "lines"),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 18),
        strip.background.x = element_rect(color = "black", 
                                          fill = c("#0099FF"),
                                          linewidth = 0),
        text = element_text(family = "Times New Roman")
      ) +
      scale_y_discrete(labels = org_fmt) +
      ggtitle(paste0(tit[i], " (n = ", format(total_count, big.mark = "", scientific = FALSE), ")")) +
      theme(plot.title = element_text(size = 18, color = "black",
                                      face = "bold", hjust = 0.5)) +
      labs(fill = "Antibiotic susceptibility testing") 
    
    plot_list[[i]] <- p
    
  }
})

# 
all_p <- wrap_plots(plot_list, nrow = 3, heights = c(0.75, 1, 1))

plot(all_p)

# Save
CairoPDF(file = "output/figure/percent_RIS_index_gpb_fungi.pdf", 
         width = 19, height = 11)
print(all_p)
dev.off()
###