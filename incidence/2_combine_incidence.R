# Clear
rm(list = ls())

suppressPackageStartupMessages({
  library(magick)
  library(pdftools)
})

# Utilities
# Safely rasterize page 1 of a single-page PDF (no Ghostscript needed)
read_pdf_as_image <- function(file_path, dpi = 300) {
  stopifnot(file.exists(file_path))
  bmp <- pdftools::pdf_render_page(file_path, page = 1, dpi = dpi)
  img <- image_read(bmp)
  # Normalize colorspace and force white background (older magick API)
  img <- image_convert(img, colorspace = "sRGB")
  img <- image_flatten(image_background(img, "white"))
  img
}

# Read, trim, add top header, then add a tiny border
process_and_label <- function(file_path, label,
                              density     = 300,
                              trim_fuzz   = 12,   # conservative to avoid over-trim
                              header_px   = 80,   # top header height
                              label_left  = 12,   # px from left within header
                              label_top   = 2,    # px from top within header
                              font_family = "Times-Bold",
                              font_size   = 45, font_weight = 800,
                              final_pad   = 6) {  # small border
  # Read page as raster
  img <- read_pdf_as_image(file_path, dpi = density)
  
  # Trim outer white margins
  img <- image_trim(img, fuzz = trim_fuzz)
  
  # Build header strip and annotate label
  info   <- image_info(img)
  header <- image_blank(width = info$width, height = header_px, color = "white")
  header <- image_annotate(
    header,
    text     = label,
    size     = font_size,
    font     = font_family,
    weight   = font_weight,
    color    = "black",
    gravity  = "northwest",
    location = sprintf("+%d+%d", label_left, label_top)
  )
  
  # Stack header over panel; ensure opaque white background
  combined <- image_append(image_join(header, img), stack = TRUE)
  combined <- image_flatten(image_background(combined, "white"))
  
  # Tiny final border
  combined <- image_border(
    combined,
    color    = "white",
    geometry = sprintf("%dx%d", final_pad, final_pad)
  )
  
  # Ensure sRGB & no alpha
  combined <- image_convert(combined, colorspace = "sRGB")
  combined <- image_flatten(image_background(combined, "white"))
  combined
}

# Bottom-pad (no scaling) to reach a target height
pad_to_height <- function(im, target_h, pad_color = "white") {
  info <- image_info(im)
  if (info$height >= target_h) return(im)
  im <- image_flatten(image_background(im, "white"))
  geom <- sprintf("%dx%d", info$width, target_h)   # width x height
  image_extent(im, geometry = geom, gravity = "northwest", color = pad_color)
}

# Make three images equal height via bottom padding and append horizontally with a white gutter
append_with_gutter <- function(img_list, gutter_px = 120) {
  stopifnot(length(img_list) == 3)
  heights  <- vapply(img_list, function(im) image_info(im)$height, numeric(1))
  target_h <- max(heights)
  
  imgs_padded <- lapply(img_list, function(im) {
    im <- image_convert(im, colorspace = "sRGB")
    im <- image_flatten(image_background(im, "white"))
    pad_to_height(im, target_h, pad_color = "white")
  })
  
  spacer <- image_blank(width = gutter_px, height = target_h, color = "white")
  frames <- list(imgs_padded[[1]], spacer, imgs_padded[[2]], spacer, imgs_padded[[3]])
  out <- image_append(do.call(image_join, frames), stack = FALSE)
  image_flatten(image_background(out, "white"))
}


# Input files (three panels per category: VAP, BSI_ICU, BSI_nonICU)
panel_files <- list(
  "All_pathogens" = c("output/incidence_meta_VAP.pdf",
                      "output/incidence_meta_BSI_ICU.pdf",
                      "output/incidence_meta_BSI_nonICU.pdf")
)


# Output configuration
out_dir <- "output/combined"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

labels     <- c("(i)", "(ii)", "(iii)")
header_px  <- 80    
gutter_px  <- 160   
label_left <- 0
label_top  <- 1

# Main loop
for (cat in names(panel_files)) {
  files <- panel_files[[cat]]
  if (length(files) != 3) {
    warning("Skipping ", cat, ": need 3 files (VAP, BSI_ICU, BSI_nonICU).")
    next
  }
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    warning("Skipping ", cat, ": missing -> ", paste(missing, collapse = ", "))
    next
  }
  
  # Process each panel: rasterize -> trim -> add header+label -> border
  img_a <- process_and_label(files[1], labels[1],
                             header_px = header_px, label_left = label_left, label_top = label_top)
  img_b <- process_and_label(files[2], labels[2],
                             header_px = header_px, label_left = label_left, label_top = label_top)
  img_c <- process_and_label(files[3], labels[3],
                             header_px = header_px, label_left = label_left, label_top = label_top)
  
  # Equalize heights via bottom padding only, then join with gutter
  triptych <- append_with_gutter(list(img_a, img_b, img_c), gutter_px = gutter_px)
  
  # Robust export: PNG first (preview), then wrap in PDF
  png_path <- file.path(out_dir, sprintf("INCIDENCE_triptych_%s.png", cat))
  image_write(triptych, path = png_path, format = "png")
  message("Wrote PNG: ", png_path)
  
  pdf_path <- file.path(out_dir, sprintf("INCIDENCE_triptych_%s.pdf", cat))
  triptych_pdf <- image_read(png_path)
  triptych_pdf <- image_convert(triptych_pdf, colorspace = "sRGB")
  image_write(triptych_pdf, path = pdf_path, format = "pdf")
  message("Wrote PDF: ", pdf_path)
}

