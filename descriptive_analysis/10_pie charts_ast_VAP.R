# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(dplyr,
                 ggplot2,
                 scatterpie,
                 RColorBrewer,
                 extrafont,
                 Cairo)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load data
df <- readRDS("data/clean data/pie_specific_anti_ast.RData")

# Load fonts
loadfonts() 

# Plot
df_plot <- df[which(df$infection_types %in% "VAP"),]
df_plot <- df_plot[which(df_plot$ris %in% c("Resistant (R)", "Intermediate (I)", "Susceptible (S)")),]
df_plot[is.na(df_plot)] <- 0

# Organisms labels
desired_order <- c(
  "Candida spp.", 
  "S. aureus",
  "Stenotrophomonas spp.",
  "Serratia/Proteus/Morganella/Enterobacter",
  "Pseudomonas spp.",
  "Acinetobacter spp.",
  "K. pneumoniae",
  "E. coli"
)

#
df_plot$org_new <- factor(df_plot$org_new, levels = desired_order)
orgi_new_labels <- levels(factor(df_plot$org_new))

# Anti groups labels
anti_group_labels <- levels(factor(df_plot$anti_group))

# Revise text
anti_group_labels[c(2, 7, 10, 13, 15)] <- 
  c("Anti-pseudomonal penicillin/\nBeta-lactamase inhibitor",
    "Fourth-generation\ncephalosporin",
    "Other Beta-lactam/\nBeta-lactamase inhibitor",
    "Sulfonamide-trimethoprim-\ncombination",
    "Third-generation\ncephalosporin")

orgi_new_labels <- sapply(orgi_new_labels, function(x) {
  if (grepl("^.*spp\\.$", x)) {
    genus <- gsub(" spp\\.$", "", x)
    return(as.expression(bquote(italic(.(genus)) ~ " spp.")))
  } else if (x == "Others") {
    return(x)
  } else {
    return(as.expression(bquote(italic(.(x)))))
  }
})

# Trans to num
df_plot$anti_group <- as.numeric(as.factor(df_plot$anti_group))
df_plot$org_new <- as.numeric(as.factor(df_plot$org_new))

#####
# Set colors
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

# 
df_plot$ris <- droplevels(df_plot$ris)
df_plot$ris <- factor(df_plot$ris ,levels = c("Resistant (R)", "Intermediate (I)", "Susceptible (S)"))
#
p <- ggplot() +
  geom_scatterpie(data = df_plot, 
                  aes(x = anti_group, y = org_new, r = 0.4, fill = anti_group),
                  cols = c(colnames(df_plot)[-(1:5)])) +
  scale_fill_manual(values = colors, guide = guide_legend(ncol = 2)) +
  scale_y_continuous(limits = c(0.5, 8.5), 
                     breaks = seq(1, 8, by = 1), 
                     expand = c(0, 0),
                     labels = orgi_new_labels) +
  scale_x_continuous(limits = c(0.5, 15.5), 
                     breaks = seq(1, 15, by = 1), 
                     expand = c(0, 0),
                     labels = anti_group_labels) +
  theme_bw() + 
  labs(x = "", y = "", title = "VAP",
       fill = "Proportion of specific antibiotic (vertical)") +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.position = "right",   
    panel.border = element_blank(),
    plot.caption = element_text(color = "black", size = 10,
                                hjust = 0, vjust = 1),
    plot.title = element_text(size = 14, color = "black",
                              face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, vjust = 0.5, 
                               colour = "black", size = 12),  
    axis.text.y = element_text(margin = margin(r = 2, t = -2),
                               colour = "black", size = 12),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    strip.text.y = element_text(size = 12, face = "bold"),
    strip.background.y = element_rect(color = "black", 
                                      fill = c("grey90"),
                                      linewidth = NA),
    text = element_text(family = "Times New Roman")
  ) + 
  facet_grid(ris ~ ., scales = "free_x", space = "free_x", switch = "x")
#
plot(p)

# Save
CairoPDF("output/figure/specific_ast_vap.pdf", width = 13.5, height = 10)
plot(p)
dev.off()
###