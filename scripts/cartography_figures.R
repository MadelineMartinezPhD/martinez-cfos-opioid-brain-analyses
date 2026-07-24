#!/usr/bin/env Rscript
# cartography_figures.R
# =============================================================================
# Figure 5, Panels A and B: cartographic node metrics mapped onto the atlas.
#   Panel A : participation coefficient (PC) coronal atlas strips + colorbar
#   Panel B : within-module degree z-score (WMDz) coronal atlas strips + colorbar
#
# Accompanies:
#   Martinez M, Ozawa A, Van Zant D, Thornberry J, Toll L.
#   Whole Brain Cellular Activation After Mu and NOP Receptor Agonism
#   Identified Differential Regional and Network Consequences.
#   [JOURNAL] 2026. [PAPER_DOI]
#
# Code author:  Madeline Martinez
# Contact:      Lawrence Toll (corresponding author)
# Affiliation:  Toll Lab, Stiles-Nicholson Brain Institute,
#               Charles E. Schmidt College of Medicine,
#               Florida Atlantic University
# Repository:   https://github.com/toll-lab-code/cfos-opioid-brain-analyses
# Archived:     https://doi.org/10.5281/zenodo.21502549
# License:      MIT (see LICENSE)
# Tested with:  R 4.4.3 (see sessionInfo.txt)
#
# Usage (from the repository root):
#   Rscript scripts/cartography_figures.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Run network_analysis.R first: this script reads the Node_Roles sheet of
# 03_Modules_and_Roles.xlsx, which that step produces.
#
# ADDITIONAL REQUIREMENTS beyond the other scripts:
#   * Scalable Brain Atlas (SBA) coronal SVG data; see the Config section and
#     the README for the expected folder contents and where to obtain it.
#   * System libraries for rsvg and magick:
#       Debian/Ubuntu  librsvg2-dev, libmagick++-dev
#       macOS/Homebrew librsvg, imagemagick
# =============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("jsonlite", "xml2", "viridisLite", "scales", "stringr",
                   "rsvg", "magick", "readxl")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(jsonlite); library(xml2); library(viridisLite); library(scales)
  library(stringr);  library(rsvg); library(magick); library(readxl)
})

# ---- Config (shared by every panel) ------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/cartography_figures.R
# Override by editing the directories below or by setting the environment
# variables CFOS_DATA_DIR / CFOS_OUTPUT_DIR before running.
DATA_ROOT   <- Sys.getenv("CFOS_DATA_DIR",   "data")
OUTPUT_ROOT <- Sys.getenv("CFOS_OUTPUT_DIR", "results")

OUTPUT_DIR <- file.path(OUTPUT_ROOT, "figures", "Figure5")
SHEET_FILE <- file.path(OUTPUT_ROOT, "03_Modules_and_Roles.xlsx")
SHEET_NAME <- "Node_Roles"
FONT       <- "Arial"
# Condition order. The atlas writes one file per condition, so this is
# order-independent here; it is kept for consistency with the other figures.
TREATMENTS <- c("Ro 64-6198", "Vehicle", "Acute Morphine", "Morphine-Dependent")

# #############################################################################
# ## SECTION A / B -- Regional atlas maps (PC, WMDz) + colorbars              ##
# #############################################################################
# =====================================================================
# Regional Atlas Visualization  (R port of generate_regional_atlas.py)
# =====================================================================
# Colors Scalable Brain Atlas (Allen Mouse Brain CCFv3, template ABA_v3) coronal
# SVG regions by per-region participation coefficient (PC) or
# within-module degree z-score (WMDz) for each treatment group.
#
# For each SBA region polygon, the metric value is resolved by:
#   1. exact acronym match in the data
#   2. else: descendants of the SBA region in the data -> averaged
#   3. else: nearest ancestor of the SBA region present in the data
#   4. else: leave the polygon gray (no coverage)
#
# Outputs HTML (interactive SVG grid) and PNG (composite raster).
#
# ---------------------------------------------------------------------
# Required packages:
#   install.packages(c("jsonlite", "xml2", "viridisLite",
#                      "scales", "stringr", "rsvg", "magick"))
# System libraries (for rsvg / magick):
#   librsvg2-dev  and  libmagick++-dev  (Debian/Ubuntu) or
#   librsvg, imagemagick (Homebrew on macOS)
# ---------------------------------------------------------------------


# -- Scalable Brain Atlas inputs ---------------------------------------
# The atlas panels need Scalable Brain Atlas coronal data for the Allen Mouse
# Brain Common Coordinate Framework v3 (SBA template ABA_v3), the same atlas
# version the cell counts are registered to. Expected in a "sba_data" folder
# containing rgb2acr.json, acr2parent.json, acr2full.json, and coronal_svg/.
# This is third-party data and is not redistributed with this repository; see
# the README for how to obtain it.
SBA_DIR <- file.path(DATA_ROOT, "sba_data")
SVG_DIR <- file.path(SBA_DIR, "coronal_svg")

# Representative coronal slices evenly spread anterior -> posterior.
# (Olfactory-bulb slices 60 & 100 dropped, plus rightmost section 500.)
SELECTED_SLICES <- c(140, 180, 220, 260, 300, 340, 380, 420, 460)

# Fraction of a slice's width that each slice overlaps its left neighbour in
# the horizontal strip (0 = no overlap, 0.12 = tuck ~12% under the previous).
OVERLAP_FRAC <- 0.38

# Color for region (and brain) outlines. "#000000" gives the black anatomical
# borders seen in published atlas figures. Set to NA to make outlines match
# each region's fill (i.e. no visible borders).
OUTLINE_COLOR <- "#000000"

# Rasterization scale relative to each slice's native SVG width. Slices are
# supersampled at this scale, then the finished strip is downscaled to the
# print size below -- so keep this comfortably above the final resolution.
RENDER_SCALE <- 3

# -- Print sizing (for assembly into a PNAS figure) --------------------
# Each strip is exported to be PLACED at this physical width and resolution.
# Layout: a 2-column PNAS figure (~7 in) with a PC column beside a WMDz column;
# ~3.3 in per strip fits two side by side with room for row labels + a gutter.
# Raise OUTPUT_DPI to 900-1200 if the strips will be enlarged beyond this width.
STRIP_WIDTH_IN  <- 3.35   # content ~82.4 mm = the A/B panel slot (side-by-side)
OUTPUT_DPI      <- 1200
STRIP_MARGIN_IN <- 0.05   # white margin around each strip (room for leaders)

# Horizontal colorbars sit under each A/B block, rendered at their FINAL
# placement size with pointsize = 6 so text is absolute pt (6 pt title, 7 pt
# ticks) and the bar drops in 1:1 -- no downscaling.
COLORBAR_WIDTH_IN  <- 2.52   # ~64 mm bar under each A/B block
COLORBAR_HEIGHT_IN <- 0.42   # 6 pt title (top) + bar + 7 pt tick labels (bottom)

METRICS <- list(
  # PC: bounded 0->1, unsigned -> sequential map (no canonical colour in the
  # lineage, where PC is encoded as node size). Swap "viridis" for "mako",
  # "magma", "rocket", etc. to taste.
  PC   = list(label = "Participation Coefficient",    cmap = "viridis", diverging = FALSE),
  # WMDz: signed z-score -> canonical blue(low)-white(0)-red(high) diverging,
  # centered at 0 (Kimbrough 2020/2021; Ardinger 2024).
  WMDz = list(label = "Within-Module Degree Z-Score", cmap = "bwr",     diverging = TRUE)
)

# -- Data source config ------------------------------------------------
# The PC / WMDz values come from the Node_Roles sheet of the network
# pipeline's output workbook (03_Modules_and_Roles.xlsx).

# Conditions to plot, in order (must match the Treatment column values).

# Which network(s) to render. One figure set is produced per sex here.
# Options present in the sheet: "Combined", "Males", "Females".
SEXES_TO_PLOT <- c("Combined", "Males", "Females")

# Color-scale groups: sexes within the same group share ONE scale per metric
# (so panels are directly comparable), and get one shared legend. Combined sits
# alone (main figure); Males+Females share a scale (supplemental sex comparison).
# Each group name is used in the legend filename. Sexes not in SEXES_TO_PLOT are
# ignored automatically.
SCALE_GROUPS <- list(
  Combined = c("Combined"),
  BySex    = c("Males", "Females")
)

# Which Node_Roles columns feed each metric.
#   Canonical / binary (matches the original labels and the Kimbrough lineage):
#       PC   -> "Participation_coef"
#       WMDz -> "WM_degree_z"
#   Weighted variants (the ones that drive the cartographic G-A roles):
#       PC   -> "Participation_coef_wt"
#       WMDz -> "WM_strength_z"
PC_COLUMN   <- "Participation_coef_wt"   # weighted PC -- matches the roles, the
WMDZ_COLUMN <- "WM_strength_z"           # network graphs, and Table S5

# -- Fixed colours -----------------------------------------------------
NO_DATA_COLOR    <- "#bfbfbf"  # regions with no metric value
BACKGROUND_COLOR <- "#ffffff"


# -- Colormap LUT (matplotlib-equivalent) ------------------------------
# matplotlib's get_cmap(name)(t) uses a 256-entry lookup table; we mirror
# that so colors match closely. t is expected in [0, 1] (vectorized).
.cmap_lut <- function(cmap_name, n = 256) {
  switch(cmap_name,
         viridis = viridisLite::viridis(n),
         magma   = viridisLite::magma(n),
         inferno = viridisLite::inferno(n),
         plasma  = viridisLite::plasma(n),
         cividis = viridisLite::cividis(n),
         mako    = viridisLite::mako(n),
         rocket  = viridisLite::rocket(n),
         # blue -> white -> red diverging (ColorBrewer RdBu, reversed so blue=low,
         # red=high). Use with a 0-centered (diverging) norm for signed metrics.
         bwr     = grDevices::colorRampPalette(c(
           "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
           "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(n),
         viridisLite::viridis(n)  # default fallback
  )
}

cmap_color <- function(cmap_name, t) {
  lut <- .cmap_lut(cmap_name)
  n <- length(lut)
  t <- pmin(pmax(t, 0), 1)            # clamp to [0,1]
  idx <- round(t * (n - 1)) + 1L
  substr(lut[idx], 1, 7)             # strip alpha if present -> #RRGGBB
}


# -- SBA atlas loading -------------------------------------------------

load_sba <- function() {
  # Read file contents explicitly, then parse. (Passing a bare path to
  # fromJSON can make it try to parse the path string itself as JSON.)
  read_json_file <- function(path, required = TRUE) {
    if (!file.exists(path)) {
      if (required) stop("Cannot find atlas file: ", path)
      warning("Optional atlas file not found; continuing without it: ", path)
      return(list())
    }
    jsonlite::fromJSON(paste(readLines(path, warn = FALSE), collapse = "\n"),
                       simplifyVector = FALSE)
  }
  
  rgb2acr <- read_json_file(file.path(SBA_DIR, "rgb2acr.json"), required = TRUE)
  names(rgb2acr) <- toupper(names(rgb2acr))
  acr2parent <- read_json_file(file.path(SBA_DIR, "acr2parent.json"), required = FALSE)
  acr2full   <- read_json_file(file.path(SBA_DIR, "acr2full.json"),   required = FALSE)
  
  # Build children index: parent -> character vector of child acronyms
  children <- list()
  for (acr in names(acr2parent)) {
    parent <- acr2parent[[acr]]
    if (is.null(parent) || is.na(parent) || !nzchar(parent)) next
    children[[parent]] <- c(children[[parent]], acr)
  }
  
  list(rgb2acr = rgb2acr, acr2parent = acr2parent,
       acr2full = acr2full, children = children)
}


all_descendants <- function(acr, children, include_self = TRUE) {
  out <- if (include_self) acr else character(0)
  stack <- children[[acr]]
  if (is.null(stack)) stack <- character(0)
  while (length(stack) > 0L) {
    a <- stack[length(stack)]
    stack <- stack[-length(stack)]
    if (a %in% out) next
    out <- c(out, a)
    kids <- children[[a]]
    if (!is.null(kids)) stack <- c(stack, kids)
  }
  out
}


ancestors <- function(acr, acr2parent) {
  out <- character(0)
  cur <- acr2parent[[acr]]
  while (!is.null(cur) && !is.na(cur) && nzchar(cur)) {
    out <- c(out, cur)
    cur <- acr2parent[[cur]]
  }
  out
}


# -- Metric value resolution -------------------------------------------

# my_values: named list/vector  acronym -> numeric (the data we have).
# Returns a closure: sba_acronym -> numeric (NA if no coverage).
# Resolution order: exact, descendants-mean, nearest-ancestor.
build_resolver <- function(my_values, acr2parent, children) {
  is_good <- vapply(my_values, function(v) !is.null(v) && !is.na(v), logical(1))
  data_keys <- names(my_values)[is_good]
  
  function(sba_acr) {
    if (sba_acr %in% data_keys) {
      return(as.numeric(my_values[[sba_acr]]))
    }
    desc <- all_descendants(sba_acr, children, include_self = FALSE)
    hit_keys <- intersect(desc, data_keys)
    if (length(hit_keys) > 0L) {
      return(mean(as.numeric(unlist(my_values[hit_keys]))))
    }
    for (anc in ancestors(sba_acr, acr2parent)) {
      if (anc %in% data_keys) {
        return(as.numeric(my_values[[anc]]))
      }
    }
    NA_real_
  }
}


# -- Color mapping -----------------------------------------------------
# In Python this returned a matplotlib Normalize object. Here we return
# a small list(vmin, vmax) plus a normalize() helper.

make_norm <- function(values, diverging) {
  arr <- as.numeric(values)
  arr <- arr[!is.na(arr)]
  if (length(arr) == 0L) return(list(vmin = 0, vmax = 1, ticks = c(0, 1)))
  if (diverging) {
    m <- max(abs(min(arr)), abs(max(arr)))
    p <- pretty(c(-m, m)); step <- if (length(p) > 1) p[2] - p[1] else m
    vmax <- ceiling(m / step) * step               # round OUT so bar ends on a tick
    return(list(vmin = -vmax, vmax = vmax, ticks = seq(-vmax, vmax, by = step)))
  }
  lo <- min(arr); hi <- max(arr)
  p <- pretty(c(lo, hi)); step <- if (length(p) > 1) p[2] - p[1] else (hi - lo)
  vmin <- floor(lo / step) * step                  # round the ends onto tick marks so the
  vmax <- ceiling(hi / step) * step                # coloured bar starts/stops at a label
  list(vmin = vmin, vmax = vmax, ticks = seq(vmin, vmax, by = step))
}

normalize_value <- function(v, norm) {
  if (norm$vmax == norm$vmin) return(rep(0.5, length(v)))
  (v - norm$vmin) / (norm$vmax - norm$vmin)
}

hex_color <- function(value, norm, cmap_name) {
  if (is.null(value) || is.na(value)) return(NO_DATA_COLOR)
  cmap_color(cmap_name, normalize_value(value, norm))
}


# -- SVG recoloring ----------------------------------------------------
# The Python version used a regex tuned to the SBA's exported format
# (one <path> per region, identical fill & stroke hex). Here we parse the
# SVG with xml2 instead: more robust to attribute ordering/whitespace,
# and we preserve the xmlns by querying with local-name() rather than
# stripping namespaces. We recolor only paths that already carry both a
# hex fill and a stroke, matching the original regex's intent.

recolor_svg <- function(svg_text, rgb2acr, resolver, norm, cmap_name) {
  doc <- read_xml(svg_text)
  
  hex6 <- "^#?[0-9A-Fa-f]{6}$"
  paths <- xml_find_all(doc, "//*[local-name()='path']")
  for (p in paths) {
    fill <- xml_attr(p, "fill")
    if (is.na(fill) || !grepl(hex6, fill)) next
    # rgb2acr keys are bare 6-char hex (no '#'); strip it before lookup.
    acr <- rgb2acr[[toupper(sub("^#", "", fill))]]
    if (is.null(acr)) next
    val <- resolver(acr)
    new_hex <- hex_color(val, norm, cmap_name)   # no data -> NO_DATA_COLOR (grey)
    xml_set_attr(p, "fill", new_hex)
    xml_set_attr(p, "stroke", if (is.na(OUTLINE_COLOR)) new_hex else OUTLINE_COLOR)
  }
  
  # Make the canvas background transparent (fill="none") so slices can be
  # overlapped without opaque rectangles occluding their neighbours.
  rects <- xml_find_all(doc, "//*[local-name()='rect']")
  for (r in rects) {
    rf <- xml_attr(r, "fill")
    if (!is.na(rf) && tolower(rf) == "#000000") {
      xml_set_attr(r, "fill", "none")
    }
  }
  
  as.character(doc)
}


# -- Rasterization -----------------------------------------------------
# cairosvg.svg2png(scale=) has no direct rsvg equivalent, so we read the
# SVG's intrinsic size and scale the requested width.

slice_to_image <- function(svg_text, scale = 1.0) {
  bytes <- charToRaw(svg_text)
  nw <- svg_native_width(svg_text)
  png_raw <- if (!is.na(nw) && nw > 0) {
    rsvg::rsvg_png(bytes, width = max(1L, round(nw * scale)))
  } else {
    rsvg::rsvg_png(bytes)  # fall back to the SVG's native size
  }
  magick::image_read(png_raw)
}

# Read the SVG's intrinsic pixel width from its width= attribute or, failing
# that, the third value of its viewBox. Replaces rsvg::rsvg_dimensions(),
# which isn't exported in all rsvg versions.
svg_native_width <- function(svg_text) {
  doc <- xml2::read_xml(svg_text)
  w <- xml2::xml_attr(doc, "width")
  if (!is.na(w)) {
    wn <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", w)))
    if (!is.na(wn) && wn > 0) return(wn)
  }
  vb <- xml2::xml_attr(doc, "viewBox")
  if (!is.na(vb)) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(vb), "[ ,]+")[[1]]))
    if (length(parts) == 4 && !is.na(parts[3]) && parts[3] > 0) return(parts[3])
  }
  NA_real_
}


# -- Composition: PNG --------------------------------------------------
# Replaces matplotlib subplot composition with magick montage.

colorbar_image <- function(norm, cmap_name, label) {
  tf <- tempfile(fileext = ".png")
  # Render at the FINAL placement size with pointsize = 6, so text is absolute:
  # cex 1 -> 6 pt (title), cex.axis 7/6 -> 7 pt (ticks). Placed 1:1 at assembly.
  grDevices::png(tf,
                 width  = round(COLORBAR_WIDTH_IN  * OUTPUT_DPI),
                 height = round(COLORBAR_HEIGHT_IN * OUTPUT_DPI),
                 res = OUTPUT_DPI, pointsize = 6)
  # Safety net: if an error occurs mid-plot, still close the device.
  ok <- FALSE
  on.exit(if (!ok && grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  
  # Horizontal bar: bold 6 pt metric title on TOP (side 3), 7 pt tick labels on
  # the BOTTOM (side 1). Sits under each A/B block in the assembly.
  op <- graphics::par(mar = c(1.6, 0.8, 1.4, 0.8))   # L/R 0.4 -> 0.8: room for end labels (0.0 / 0.7)
  n <- 256L
  cols <- cmap_color(cmap_name, (0:(n - 1)) / (n - 1))
  graphics::image(
    x = seq(norm$vmin, norm$vmax, length.out = n), y = 1,
    z = matrix(seq_len(n), nrow = n, ncol = 1),
    col = cols, axes = FALSE, xlab = "", ylab = ""
  )
  graphics::axis(1, at = norm$ticks, cex.axis = 7/6, mgp = c(3, 0.4, 0), tcl = -0.3)
  graphics::mtext(label, side = 3, line = 0.3, font = 2, cex = 1)
  graphics::par(op)
  
  grDevices::dev.off()   # close & flush to disk BEFORE reading it back
  ok <- TRUE
  magick::image_read(tf)
}

compose_png <- function(slice_imgs, png_path,
                        overlap_frac = OVERLAP_FRAC, bg = "white") {
  # Trim each slice's transparent margin so neighbours tuck together by their
  # actual brain content rather than their bounding boxes.
  trimmed <- lapply(slice_imgs, magick::image_trim)
  info <- lapply(trimmed, magick::image_info)
  ws <- vapply(info, function(x) x$width,  numeric(1))
  hs <- vapply(info, function(x) x$height, numeric(1))
  H  <- max(hs)
  
  # Horizontal overlap, sized from the mean slice width.
  ov <- round(overlap_frac * mean(ws))
  
  # Left edge x for each slice (each starts ov px before its predecessor ends).
  xs <- numeric(length(trimmed))
  cur <- 0
  for (i in seq_along(trimmed)) {
    xs[i] <- cur
    cur <- cur + ws[i] - ov
  }
  total_w <- cur + ov
  
  # Composite onto a solid background canvas. Slices keep transparent margins,
  # so overlaps reveal the neighbour beneath rather than an opaque box; later
  # (posterior) slices land on top.
  strip <- magick::image_blank(total_w, H, color = bg)
  for (i in seq_along(trimmed)) {
    yoff <- round((H - hs[i]) / 2)
    strip <- magick::image_composite(
      strip, trimmed[[i]], offset = sprintf("+%d+%d", round(xs[i]), yoff)
    )
  }
  
  # Downscale the supersampled strip to the target print width, add the margin,
  # and tag the PNG resolution so it imports at STRIP_WIDTH_IN @ OUTPUT_DPI.
  margin_px <- round(STRIP_MARGIN_IN * OUTPUT_DPI)
  inner_px  <- max(1L, round(STRIP_WIDTH_IN * OUTPUT_DPI) - 2L * margin_px)
  strip <- magick::image_resize(strip, as.character(inner_px))
  if (margin_px > 0) {
    strip <- magick::image_border(strip, bg, sprintf("%dx%d", margin_px, margin_px))
  }
  magick::image_write(strip, path = png_path, format = "png",
                      density = as.character(OUTPUT_DPI))
}

# Write a standalone horizontal colorbar (legend) for one metric's shared scale,
# placed 1:1 under its A/B strip block at assembly.
write_legend <- function(norm, cmap_name, label, png_path) {
  cb <- colorbar_image(norm, cmap_name, label)
  magick::image_write(cb, path = png_path, format = "png",
                      density = as.character(OUTPUT_DPI))
}


# -- Data loader (reads the Node_Roles sheet) --------------------------
# Builds the nested structure the main loop expects:
#   study(= Sex) -> group(= Treatment) -> list(
#       participation_coefficient = named list  region -> PC,
#       wmdz                      = named list  module_id -> (region -> WMDz)
#   )
# WMDz is kept nested by Network_module to mirror the original Python
# contract; flatten_wmdz() collapses it to region -> value downstream.

load_all_studies <- function() {
  if (!file.exists(SHEET_FILE)) {
    stop("Cannot find data workbook: ", SHEET_FILE)
  }
  roles <- as.data.frame(read_excel(SHEET_FILE, sheet = SHEET_NAME),
                         stringsAsFactors = FALSE)
  
  needed <- c("Treatment", "Sex", "Region", "Network_module",
              PC_COLUMN, WMDZ_COLUMN)
  missing <- setdiff(needed, names(roles))
  if (length(missing) > 0) {
    stop("Sheet '", SHEET_NAME, "' is missing column(s): ",
         paste(missing, collapse = ", "))
  }
  
  roles$Region <- trimws(as.character(roles$Region))
  roles <- roles[!is.na(roles$Region) & nzchar(roles$Region), ]
  
  studies <- list()
  
  for (sx in SEXES_TO_PLOT) {
    sub_sex <- roles[roles$Sex == sx, , drop = FALSE]
    if (nrow(sub_sex) == 0L) {
      message("  Warning: no rows for Sex == '", sx, "'")
      next
    }
    groups <- list()
    for (trt in TREATMENTS) {
      sub <- sub_sex[sub_sex$Treatment == trt, , drop = FALSE]
      if (nrow(sub) == 0L) {
        message("  Warning: no rows for Treatment == '", trt, "' (Sex ", sx, ")")
        next
      }
      
      # participation_coefficient: region -> value
      pc <- as.list(suppressWarnings(as.numeric(sub[[PC_COLUMN]])))
      names(pc) <- sub$Region
      
      # wmdz: module_id -> (region -> value)
      wmdz <- list()
      wm_vals <- suppressWarnings(as.numeric(sub[[WMDZ_COLUMN]]))
      for (i in seq_len(nrow(sub))) {
        mod <- as.character(sub$Network_module[i])
        if (is.null(wmdz[[mod]])) wmdz[[mod]] <- list()
        wmdz[[mod]][[sub$Region[i]]] <- wm_vals[i]
      }
      
      groups[[trt]] <- list(participation_coefficient = pc, wmdz = wmdz)
    }
    if (length(groups) > 0L) studies[[sx]] <- groups
  }
  
  if (length(studies) == 0L) {
    stop("No data assembled -- check SEXES_TO_PLOT / TREATMENTS against the sheet.")
  }
  studies
}

safe_filename <- function(s) {
  s <- gsub("[^A-Za-z0-9._-]+", "_", s)
  gsub("^_+|_+$", "", s)
}


# -- Main --------------------------------------------------------------

# Each region belongs to one cluster; merge to {region: value}.
flatten_wmdz <- function(wmdz_by_cluster) {
  out <- list()
  for (cid in names(wmdz_by_cluster)) {
    region_map <- wmdz_by_cluster[[cid]]
    for (region in names(region_map)) {
      out[[region]] <- region_map[[region]]
    }
  }
  out
}

main <- function() {
  sba <- load_sba()
  rgb2acr    <- sba$rgb2acr
  acr2parent <- sba$acr2parent
  acr2full   <- sba$acr2full
  children   <- sba$children
  
  # Cache SVGs for selected slices once.
  slice_svgs <- list()
  for (idx in SELECTED_SLICES) {
    path <- file.path(SVG_DIR, sprintf("Annotation2014_141_%04d.svg", idx))
    slice_svgs[[as.character(idx)]] <-
      paste(readLines(path, warn = FALSE), collapse = "\n")
  }
  
  all_data <- load_all_studies()
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # SBA acronyms that actually appear in the selected slices. This depends
  # only on the SVGs, so compute it once and reuse for every group/metric.
  acrs_in_view <- character(0)
  for (svg in slice_svgs) {
    fills <- str_match_all(svg, 'fill="#([0-9A-Fa-f]{6})"')[[1]]
    if (nrow(fills) > 0L) {
      for (f in fills[, 2]) {
        acr <- rgb2acr[[toupper(f)]]
        if (!is.null(acr)) acrs_in_view <- c(acrs_in_view, acr)
      }
    }
  }
  acrs_in_view <- unique(acrs_in_view)
  
  # Determine the single-vs-multi-sex suffix policy: when only one sex is being
  # plotted, filenames stay clean (no sex tag); otherwise every strip is tagged.
  sexes_present <- names(all_data)
  tag_sex <- length(sexes_present) > 1L
  
  for (grp_label in names(SCALE_GROUPS)) {
    sexes <- intersect(SCALE_GROUPS[[grp_label]], sexes_present)
    if (length(sexes) == 0L) next
    
    for (metric_key in names(METRICS)) {
      meta <- METRICS[[metric_key]]
      
      # Build a resolver for every (sex x condition) in this scale group, and
      # pool ALL their values so the colour scale is shared across the group.
      resolvers   <- list()        # keyed "sex||condition"
      pooled_vals <- numeric(0)
      for (sx in sexes) {
        for (group_name in names(all_data[[sx]])) {
          gd <- all_data[[sx]][[group_name]]
          my_values <- if (metric_key == "PC") gd$participation_coefficient
          else flatten_wmdz(gd$wmdz)
          resolver <- build_resolver(my_values, acr2parent, children)
          resolvers[[paste(sx, group_name, sep = "||")]] <- resolver
          for (acr in acrs_in_view) {
            v <- resolver(acr)
            if (!is.na(v)) pooled_vals <- c(pooled_vals, v)
          }
        }
      }
      
      # ONE norm for the whole scale group, applied to every panel in it.
      norm <- make_norm(pooled_vals, meta$diverging)
      cat(sprintf("[%s / %s] shared scale across {%s}: [%.3f, %.3f]\n",
                  grp_label, metric_key, paste(sexes, collapse = ", "),
                  norm$vmin, norm$vmax))
      
      # One shared legend per scale group + metric.
      legend_base <- sprintf("regional_%s_%s_legend", metric_key, grp_label)
      write_legend(norm, meta$cmap, meta$label,
                   file.path(OUTPUT_DIR, paste0(legend_base, ".png")))
      cat(sprintf("  wrote %s.png  (shared legend)\n", legend_base))
      
      # Render each (sex, condition) strip with the group's shared scale.
      for (sx in sexes) {
        for (group_name in names(all_data[[sx]])) {
          resolver <- resolvers[[paste(sx, group_name, sep = "||")]]
          
          slice_imgs <- list()
          for (idx in SELECTED_SLICES) {
            svg_new <- recolor_svg(
              slice_svgs[[as.character(idx)]], rgb2acr, resolver, norm, meta$cmap
            )
            slice_imgs[[length(slice_imgs) + 1L]] <- slice_to_image(svg_new, scale = RENDER_SCALE)
          }
          n_colored <- sum(vapply(acrs_in_view,
                                  function(a) !is.na(resolver(a)), logical(1)))
          
          base <- if (tag_sex)
            sprintf("regional_%s_%s_%s", metric_key, safe_filename(group_name),
                    safe_filename(sx))
          else
            sprintf("regional_%s_%s", metric_key, safe_filename(group_name))
          compose_png(slice_imgs, file.path(OUTPUT_DIR, paste0(base, ".png")))
          
          cat(sprintf("  wrote %s.png  (n=%d colored regions)\n",
                      base, n_colored))
        }
      }
    }
  }
  
  cat("Done.\n")
}

# Generate the atlas maps for Panels A and B.
main()
