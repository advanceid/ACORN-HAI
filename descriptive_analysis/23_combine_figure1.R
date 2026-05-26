# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(magick, pdftools, png, Cairo, grid)
})

# Read page 1 of PDF as image
read_pdf_as_image <- function(file_path, dpi = 600) {
  stopifnot(file.exists(file_path))
  
  bmp <- pdftools::pdf_render_page(
    pdf  = file_path,
    page = 1,
    dpi  = dpi
  )
  
  img <- image_read(bmp)
  img <- image_convert(img, colorspace = "sRGB")
  img <- image_flatten(image_background(img, "white"))
  img
}

# Trim outer white margins
trim_pdf_image <- function(img, trim_fuzz = 8) {
  img <- image_trim(img, fuzz = trim_fuzz)
  img <- image_convert(img, colorspace = "sRGB")
  img <- image_flatten(image_background(img, "white"))
  img
}


# Add label AFTER resize
add_panel_label <- function(img, label,
                            top_pad     = 140,
                            left_pad    = 70,
                            right_pad   = 20,
                            bottom_pad  = 20,
                            label_x     = 18,
                            label_y     = 12,
                            font_family = "Times",
                            font_size   = 90,
                            font_weight = 700) {
  
  info <- image_info(img)
  
  canvas_w <- info$width + left_pad + right_pad
  canvas_h <- info$height + top_pad + bottom_pad
  
  # Put original image at bottom-right so top/left white space is reserved for label
  out <- image_extent(
    img,
    geometry = sprintf("%dx%d", canvas_w, canvas_h),
    gravity  = "southeast",
    color    = "white"
  )
  
  out <- image_annotate(
    out,
    text     = label,
    size     = font_size,
    font     = font_family,
    weight   = font_weight,
    color    = "black",
    gravity  = "northwest",
    location = sprintf("+%d+%d", label_x, label_y)
  )
  
  out <- image_convert(out, colorspace = "sRGB")
  out <- image_flatten(image_background(out, "white"))
  out
}


# Pad to target width for vertical stacking
pad_to_width <- function(im, target_w, pad_color = "white") {
  info <- image_info(im)
  
  if (info$width >= target_w) return(im)
  
  image_extent(
    im,
    geometry = sprintf("%dx%d", target_w, info$height),
    gravity  = "northwest",
    color    = pad_color
  )
}


# Combine 2 images vertically with white gutter
append_vertical_with_gutter <- function(img_list, gutter_px = 80) {
  stopifnot(length(img_list) == 2)
  
  widths   <- vapply(img_list, function(im) image_info(im)$width, numeric(1))
  target_w <- max(widths)
  
  imgs_pad <- lapply(img_list, function(im) {
    im <- image_convert(im, colorspace = "sRGB")
    im <- image_flatten(image_background(im, "white"))
    pad_to_width(im, target_w, pad_color = "white")
  })
  
  spacer <- image_blank(
    width  = target_w,
    height = gutter_px,
    color  = "white"
  )
  
  out <- image_append(
    do.call(image_join, list(imgs_pad[[1]], spacer, imgs_pad[[2]])),
    stack = TRUE
  )
  
  out <- image_convert(out, colorspace = "sRGB")
  out <- image_flatten(image_background(out, "white"))
  out
}


# Input files
file1 <- "output/figure/pathogen_age_proportion.pdf"
file2 <- "output/figure/proportion_resist_top3_h.pdf"

# Read original PDFs
img_a <- read_pdf_as_image(file1, dpi = 600)
img_b <- read_pdf_as_image(file2, dpi = 600)

# Trim first
img_a <- trim_pdf_image(img_a, trim_fuzz = 8)
img_b <- trim_pdf_image(img_b, trim_fuzz = 4)

# Resize figure 2 first, BEFORE adding label
# Match width of figure 1
w1 <- image_info(img_a)$width
img_b <- image_resize(img_b, paste0(w1, "x"))

# Now add labels AFTER resize
img_a <- add_panel_label(
  img_a,
  label       = "(a)",
  top_pad     = 140,
  left_pad    = 70,
  right_pad   = 20,
  bottom_pad  = 20,
  label_x     = 18,
  label_y     = 12,
  font_family = "Times",
  font_size   = 90,
  font_weight =1000
)

img_b <- add_panel_label(
  img_b,
  label       = "(b)",
  top_pad     = 80,
  left_pad    = 70,
  right_pad   = 20,
  bottom_pad  = 20,
  label_x     = 18,
  label_y     = 12,
  font_family = "Times",
  font_size   = 90,
  font_weight = 1000
)


# Combine vertically
combined_img <- append_vertical_with_gutter(
  img_list  = list(img_a, img_b),
  gutter_px = 80
)

# Make sure output folder exists
dir.create("output/figure", recursive = TRUE, showWarnings = FALSE)

# Export PNG first
png_path <- file.path("output/figure/com_figure1_vertical.png")
pdf_path <- file.path("output/figure/com_figure1_vertical.pdf")

combined_img <- image_convert(combined_img, colorspace = "sRGB")
combined_img <- image_flatten(image_background(combined_img, "white"))

image_write(
  combined_img,
  path = png_path,
  format = "png"
)

# Read PNG back
img_png <- png::readPNG(png_path)

# Get image size
info <- image_info(combined_img)
w_in <- info$width / 600
h_in <- info$height / 600

# Export to PDF using Cairo
CairoPDF(
  file = pdf_path,
  width = w_in,
  height = h_in
)

grid::grid.newpage()
grid::grid.raster(
  img_png,
  x = 0,
  y = 1,
  width = 1,
  height = 1,
  just = c("left", "top"),
  interpolate = TRUE
)

dev.off()

#####
