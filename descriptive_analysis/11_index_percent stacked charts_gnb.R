# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magrittr, 
                 dplyr,
                 ggplot2,
                 grid,
                 gridExtra,
                 patchwork,
                 extrafont,
                 Cairo,
                 ggtext)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean_data/percent_RIS.RData")
index_number <- readRDS("data/clean_data/organism_counts.RData")

# Load fonts
loadfonts()

# Make plots
data_list <- plot_list <- list()

# Split and clean data
system.time({
  for (i in seq_along(df)) {
    
    df_used <- df[[i]][which(df[[i]]$organism_group == "GNB"),]
    df_used$organism_names <- droplevels(df_used$organism_names)
    df_used$anti_group <- droplevels(df_used$anti_group)
    
    # New order
    new_order <- c("Fluoroquinolone",
                   "Aminoglycoside",                                        
                   "Carbapenem", 
                   "Fourth-generation cephalosporin",                      
                   "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
                   "Other beta-lactam/beta-lactamase inhibitor",
                   "Third-generation cephalosporin",
                   "Sulfonamide-trimethoprim-combination",
                   "Tetracycline")
    
    df_used$anti_group <- factor(df_used$anti_group, levels = new_order)
    
    # Rename levels
    levels(df_used$anti_group)[4:8] <- c(
      "Fourth-generation<br>cephalosporin",
      "Anti-pseudomonal<br>penicillin/beta-<br>lactamase inhibitor<sup>*</sup>",
      "Other beta-lactam/<br>beta-lactamase<br>inhibitor<sup>#</sup>",
      "Third-generation<br>cephalosporin",
      "Sulfonamide-<br>trimethoprim-<br>combination"
    )
    
    data_list[[i]] <- df_used
  }
})

# Swap 1st and 3rd plots
temp <- data_list[[1]]
data_list[[1]] <- data_list[[3]]
data_list[[3]] <- temp

# Unify anti_group levels across all datasets
all_levels <- levels(data_list[[1]]$anti_group)
data_list <- lapply(data_list, function(df) {
  df$anti_group <- factor(df$anti_group, levels = all_levels)
  df
})

# Plot
options(expressions = 50000)

system.time({
  for (i in 1:3) {
    df_used <- data_list[[i]]
    
    # Generate display labels
    org_labels <- as.character(df_used$organism_names)
    org_fmt <- sapply(org_labels, function(x) {
      if (x == "Serratia/Proteus/Morganella/Enterobacter") {
        return(bquote(atop(italic("Serratia/Proteus/"), 
                           italic("Morganella/Enterobacter"))))
      } else if (x %in% c("Others")) {
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
    
    tit <- levels(df_used$infection_types)
    total_count <- index_number$GNB[i]
    df_used$organism_names <- factor(df_used$organism_names, 
                                     levels = unique(df_used$organism_names))
    
    p <- ggplot(df_used, 
                aes(y = organism_names,
                    x = per_count, 
                    fill = anti_susceptibility)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      geom_text(aes(label = sum_count, x = 0.9, y = organism_names), 
                vjust = 0.5, size = 4.5, fontface = "bold", color = "white") +
      facet_wrap(~ anti_group, nrow = 1, drop = FALSE) +
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
        strip.text.x = ggtext::element_markdown(size = 16),
        strip.background.x = element_rect(color = "black", 
                                          fill = c("#0099FF"),
                                          linewidth = 0),
        text = element_text(family = "Times New Roman")
      ) + 
      scale_y_discrete(labels = org_fmt) +
      ggtitle(paste0(tit[i], " (n = ", format(total_count, big.mark = "", scientific = FALSE), ")")) +
      theme(plot.title = element_text(size = 18, color = "black",
                                      face = "bold", hjust = 0.5)) +
      guides(fill = guide_legend(title = "Antibiotic susceptibility testing", nrow = 1))
    
    plot_list[[i]] <- p
  }
})

# Combine and plot
all_p <- wrap_plots(plot_list, nrow = 3)
plot(all_p)

# Save
CairoPDF(file = "output/figure/percent_RIS_index_gnb_ppt.pdf", 
         width = 23, height = 15)
print(all_p)
dev.off()
###
