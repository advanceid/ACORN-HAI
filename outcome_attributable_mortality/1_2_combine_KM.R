# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magick)
})

# Set working directory
wd <- "./"
setwd(wd)

# Create file names
pdf_files <- paste0("output/sub_figure/att28_death_KM_0", 1:6, "_overall.pdf")

# Read and convert first page of each PDF to image
pdf_imgs <- lapply(pdf_files, function(path) image_read_pdf(path, pages = 1))

# Ensure all are single-page images (flatten)
pdf_imgs <- lapply(pdf_imgs, function(img) image_flatten(img))

# Combine row-wise: 3 columns × 2 rows
row1 <- image_append(c(pdf_imgs[[1]], pdf_imgs[[2]], pdf_imgs[[3]]), stack = FALSE)
row2 <- image_append(c(pdf_imgs[[4]], pdf_imgs[[5]], pdf_imgs[[6]]), stack = FALSE)
final_img <- image_append(c(row1, row2), stack = TRUE)

# Save
image_write(final_img, path = "output/att28_death_KM_combined_overall.pdf", format = "pdf")
###
