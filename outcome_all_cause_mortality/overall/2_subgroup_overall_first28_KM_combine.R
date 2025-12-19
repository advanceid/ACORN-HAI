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

# Create a vector of PDF file names
pdf_files <- c("output/pdf/first28_death_KM_overall.pdf",
               "output/pdf/first28_death_KM.pdf")

# Read the PDFs into a list of images
pdf_images <- lapply(pdf_files, image_read_pdf)

# Define the labels
labels <- c("(a) Overall",
            "(b) Subgroups")

# Annotate each image with corresponding label
pdf_images_annotated <- mapply(function(image, label) {
  image_annotate(image, label, gravity = "northwest", size = 10, color = "black", location = "+1+20", weight = 400, font = "Times New Roman")
}, pdf_images, labels, SIMPLIFY = FALSE)

# Combine 
final_image <- image_append(c(pdf_images_annotated[[1]], 
                              pdf_images_annotated[[2]]), stack = T)

# 
print(final_image)

# Save
image_write(final_image, path = "output/pdf/first28_death_KM_com.pdf", format = "pdf")
###