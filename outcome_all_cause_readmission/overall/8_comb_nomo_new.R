# Clear environment
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magick, grid, ragg)
})

# Set working directory
setwd("./")

# Load fonts
loadfonts()

# Step 1: Read and label two PDF images
pdf1_image <- image_read_pdf("output/pdf/nomo_readmission.pdf")
pdf2_image <- image_read_pdf("output/pdf/nomo_death.pdf")

pdf1_image <- image_annotate(pdf1_image, "(a)", gravity = "northwest", size = 18,
                             color = "black", location = "+220+40", weight = 400,
                             font = "Times New Roman")
pdf2_image <- image_annotate(pdf2_image, "(b)", gravity = "northwest", size = 18,
                             color = "black", location = "+220+40", weight = 400,
                             font = "Times New Roman")

# Combine vertically
main_image <- image_append(c(pdf1_image, pdf2_image), stack = TRUE)

# Step 2: Dynamically get width of main image
main_info <- image_info(main_image)
note_width <- main_info$width
note_height <- 250 # Give enough height for 2 lines

# Step 3: Create annotation image (2 lines, with italic)
note_path <- "output/pdf/note_image.png"
ragg::agg_png(note_path, width = note_width, height = note_height, res = 300, bg = "white")

grid.text(
  expression("Abbreviations: INF = Infectious Disease, GIT = Gastrointestinal Disorder, PMD = Pulmonary Disease, OTH = Others, VAP = Ventilator-associated Pneumonia, BSI = Bloodstream Infection, CRA = Carbapenem-resistant"),
  x = 0.03, y = 0.85, just = "left",
  gp = gpar(fontsize = 18, fontfamily = "Times New Roman")
)

grid.text(
  expression(italic(Acinetobacter)~" spp.; 3GCRE = Third-generation Cephalosporins-resistant Enterobacterales, CRE = Carbapenem-resistant Enterobacterales, ICU = Intensive Care Unit, HD = High Dependency."),
  x = 0.097, y = 0.5, just = "left",
  gp = gpar(fontsize = 18, fontfamily = "Times New Roman")
)

dev.off()

# Step 4: Read note image and combine
note_image <- image_read(note_path)
final_image <- image_append(c(main_image, note_image), stack = TRUE)

# Step 5: Save final image as PDF
image_write(final_image, path = "output/pdf/com_nomo.pdf", format = "pdf")
