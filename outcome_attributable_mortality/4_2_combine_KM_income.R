# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magick)
})

# Define working directory
wd <- "./"
setwd(wd)

# Create a vector of PDF file names
pdf_files <- paste0("output/sub_figure/att28_death_KM_0", 1:6, "_income.pdf")

# Read the PDFs into a list of images
pdf_images <- lapply(pdf_files, image_read_pdf)

# Define the labels
labels <- c("(a)",
            "(b)",
            "(c)",
            "(d)",
            "(e)",
            "(f)")

# Annotate each image with corresponding label
pdf_images_annotated <- mapply(function(image, label) {
  image_annotate(image, label, gravity = "northwest", size = 18, color = "black", location = "+20+1", weight = 400, font = "Times New Roman")
}, pdf_images, labels, SIMPLIFY = FALSE)

# Combine 1–3
final_image_1 <- image_append(c(pdf_images_annotated[[1]],
                                pdf_images_annotated[[2]],
                                pdf_images_annotated[[3]]), stack = TRUE)

# Combine 4–6
final_image_2 <- image_append(c(pdf_images_annotated[[4]],
                                pdf_images_annotated[[5]],
                                pdf_images_annotated[[6]]), stack = TRUE)

# Save
image_write(final_image_1, path = "output/com_att28_KM_1to3_income.pdf", format = "pdf")
image_write(final_image_2, path = "output/com_att28_KM_4to6_income.pdf", format = "pdf")

