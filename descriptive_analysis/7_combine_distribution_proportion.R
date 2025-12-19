# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magick,
                 extrafont)
})

# Define working directory
wd <- "./"
setwd(wd)

# Load fonts
loadfonts()

# pdf to fig
pdf1_image <- image_read_pdf("output/figure/proportion_index_episode.pdf")
pdf2_image <- image_read_pdf("output/figure/index_organism_distribution_v.pdf")
pdf3_image <- image_read_pdf("output/figure/index_pathogen_types_distribution_v.pdf")

# Combine
combined_image <- image_append(c(pdf1_image[1], pdf3_image[1], pdf2_image[1]))

# Save final image
image_write(combined_image, path = "output/figure/com_proportion_distribution_v.pdf", format = "pdf")
###