# regional_figures.R
# =============================================================================
# Figures 2 and 3: regional activation results.
#
# Accompanies:
#   Martinez M, Ozawa A, Van Zandt D, Thornberry J, Toll L.
#   Whole Brain Cellular Activation After Mu and NOP Receptor Agonism
#   Identified Differential Regional and Network Consequences.
#   [JOURNAL] 2026. [PAPER_DOI]
#
# Code author:  Madeline Martinez
# Contact:      Lawrence Toll (corresponding author)
# Affiliation:  Toll Lab, Stiles-Nicholson Brain Institute,
#               Charles E. Schmidt College of Medicine,
#               Florida Atlantic University
# Repository:   https://github.com/MadelineMartinezPhD/martinez-cfos-opioid-brain-analyses
# Archived:     [ZENODO_CODE_DOI]
# License:      MIT (see LICENSE)
# Tested with:  R 4.4.3 (see sessionInfo.txt)
#
# Usage (from the repository root):
#   Rscript scripts/regional_figures.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Run regional_analysis.R first: this script reads the statistics workbooks that
# it produces. Panels are rendered individually and assembled into the final
# figures separately.
#
# PANEL MANIFEST
# Figure 2, Panel A: Lollipop - Veh vs Acute Morphine; top 30 regions by Cohen's
#                    |d| drawn from the 74 significant at uncorrected p<0.05.
#                    FDR-significant (q<0.05) regions = filled points; uncorrected-
#                    only = open points. Table S2.
# Figure 2, Panel B: Significant-region count summary - uncorrected & FDR counts
#                    per drug condition (combined sex)
# Figure 2, Panel C: Whole-brain total cFos+ cells per animal, by condition
#                    (combined sex; one-way ANOVA)
# Figure 2, Panel D: Representative light-sheet images - assembled in BioRender,
#                    not produced by this script
# Figure 2, Panel E: Faceted exemplar bars - MBO, ASO, TRS, MEA; all 4 conditions;
#                    combined sex
# Figure 3, Panel A: Diverging bar chart - males vs. females, Veh vs. Acute Morphine
# Figure 3, Panel B: Faceted sex-stratified bars - SNr/AVPV/LC (rows) x
#                    Males/Females (cols); shared y per region
# Figure 3, Panel C: Acute vs. dependent effect-size comparison - Veh vs Acute |d|
#                    vs Veh vs Dependent |d|; descriptive FDR-vs-vehicle grouping
# Figure 3, Panel D: Faceted acute-vs-dependent exemplar bars - DMH (dependence-
#                    only), NOD (acute-only); combined sex, all 4 conditions
#
# Output: 1200 dpi LZW TIFF + matching PNG per panel, for BioRender assembly
# =============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("tidyverse", "readxl", "showtext", "ggbeeswarm", "scales")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(tidyverse)
library(readxl)
library(showtext)
library(ggbeeswarm)
library(scales)

# ---- Font --------------------------------------------------------------------
# Figures use Arial. The default search path is the macOS system font directory;
# set the ARIAL_DIR environment variable to point somewhere else. If the font
# files are not found, the script falls back to the device default sans family,
# which alters glyph metrics slightly but nothing else about the figures.
arial_dir   <- Sys.getenv("ARIAL_DIR", "/System/Library/Fonts/Supplemental")
arial_files <- file.path(arial_dir, c("Arial.ttf", "Arial Bold.ttf",
                                      "Arial Italic.ttf", "Arial Bold Italic.ttf"))
FONT_FAMILY <- "sans"
if (all(file.exists(arial_files))) {
  font_add("Arial", regular = arial_files[1], bold = arial_files[2],
           italic = arial_files[3], bolditalic = arial_files[4])
  FONT_FAMILY <- "Arial"
} else {
  message("Arial not found in ", arial_dir, ". Falling back to the default sans ",
          "family; set ARIAL_DIR to the directory containing Arial.ttf to ",
          "reproduce the published figures exactly.")
}
showtext_auto()
showtext_opts(dpi = 1200)

# ---- Config ------------------------------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/regional_figures.R
# Override by editing the directories below or by setting the environment
# variables CFOS_DATA_DIR / CFOS_OUTPUT_DIR / CFOS_SCRIPT_DIR before running.
DATA_DIR   <- Sys.getenv("CFOS_DATA_DIR",   "data")
OUTPUT_DIR <- Sys.getenv("CFOS_OUTPUT_DIR", "results")
SCRIPT_DIR <- Sys.getenv("CFOS_SCRIPT_DIR", "scripts")

raw_file   <- file.path(DATA_DIR, "Drug_Composite.xlsx")
vol_file   <- file.path(DATA_DIR, "ccfv3_volumes.xlsx")

# Figure 2 panels save here; reassigned to Figure3/ before the Figure 3 section
output_dir <- file.path(OUTPUT_DIR, "figures", "Figure2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -- Colors --------------------------------------------------------------------
# The division palette and the region-to-division lookup are both defined in
# divisions.R, so the 13-division scheme is identical across every figure.
source(file.path(SCRIPT_DIR, "divisions.R"))
structure_colors <- division_colors            # 13 divisions, canonical palette

# Canonical region -> division lookup (Allen CCFv3 divisions, with the lab's two
# functional groupings: amygdalar nuclei kept as one Amygdala division; lateral
# septum grouped with Pallidum). Every panel resolves `structure` from THIS table
# by acronym, so the grouping cannot drift between figures. Generated from the
# Node_Roles Structure column (region_division_lookup.csv).
division_lookup_file <- file.path(DATA_DIR, "region_division_lookup.csv")
.divtab <- read.csv(division_lookup_file, stringsAsFactors = FALSE)
div_map <- setNames(.divtab$Structure, .divtab$Region)

condition_colors <- c(
  "Vehicle"            = "#9e9e9e",   # grey (control)
  "Acute Morphine"     = "#c98ba0",   # dusty rose
  "Morphine-Dependent" = "#8f72b3",   # dusty mauve-purple
  "Ro 64-6198"         = "#4fae9f"    # soft teal-green
)

condition_order <- c("Vehicle", "Acute Morphine", "Morphine-Dependent", "Ro 64-6198")

sex_colors <- c("Males" = "#3f6fa3", "Females" = "#7d2452")   # muted blue / wine


# -- Save helper -------------------------------------------------------------
# Writes each panel twice at 1200 dpi: a TIFF (LZW) for archival/submission and
# a PNG for BioRender import. PNG is lossless and carries the dpi, so panels land
# at their true physical size in BioRender (no JPEG artifacts, no mis-scaling).
# Pass the .tiff path; the .png path is derived from it.
save_panel <- function(plot, path_tiff, width, height) {
  ggsave(path_tiff, plot, width = width, height = height, units = "mm",
         dpi = 1200, compression = "lzw", device = "tiff")
  ggsave(sub("\\.tiff?$", ".png", path_tiff), plot,
         width = width, height = height, units = "mm",
         dpi = 1200, device = "png")
}


# =============================================================================
# PANEL A - Lollipop: top 30 regions by Cohen's |d| among those significant
#           at uncorrected p<0.05 (Vehicle vs. Acute Morphine, combined sex).
#           Points are filled where the region also survives FDR (q<0.05) and
#           open where it is significant at uncorrected p<0.05 only (e.g. AVPV,
#           q=0.06). Full list of all 74 significant regions is in Table S2.
# =============================================================================

fig2a_data <- tibble(
  acronym = c("MBO","ASO","PVpo","TRS","NOD","SBPV","NDB","AVPV","LPO","AHN",
              "ATN","PVa","LHA","RAmb","PVi","NLL","MEA","SNr","MPN","PGRN",
              "PVH","PeF","RT","SI","AT","TT","VeCB","GPe","PRM","PN"),
  label = c(
    "Mammillary body",
    "Accessory supraoptic group",
    "PVN, preoptic part",
    "Triangular nucleus of septum",
    "Nodulus",
    "Subparaventricular zone",
    "Diagonal band nucleus",
    "Anteroventral periventricular nucleus",
    "Lateral preoptic area",
    "Anterior hypothalamic nucleus",
    "Anterior group, dorsal thalamus",
    "PVN, anterior part",
    "Lateral hypothalamic area",
    "Midbrain raphe nuclei",
    "PVN, intermediate part",
    "Nucleus of the lateral lemniscus",
    "Medial amygdalar nucleus",
    "Substantia nigra, reticular part",
    "Medial preoptic nucleus",
    "Paragigantocellular reticular nucleus",
    "Paraventricular hypothalamic nucleus",
    "Perifornical nucleus",
    "Reticular nucleus of thalamus",
    "Substantia innominata",
    "Anterior tegmental nucleus",
    "Taenia tecta",
    "Vestibulocerebellar nucleus",
    "Globus pallidus, external",
    "Paramedian lobule",
    "Paranigral nucleus"
  ),
  abs_d = c(2.070,1.945,1.775,1.747,1.742,1.706,1.702,1.663,1.545,1.530,
            1.523,1.523,1.504,1.464,1.459,1.439,1.431,1.419,1.418,1.398,
            1.389,1.380,1.360,1.357,1.325,1.309,1.304,1.294,1.288,1.281),
  structure = c(
    "Hypothalamus","Hypothalamus","Hypothalamus","Pallidum/Septum",
    "Hindbrain","Hypothalamus","Pallidum/Septum","Hypothalamus",
    "Hypothalamus","Hypothalamus","Thalamus","Hypothalamus",
    "Hypothalamus","Midbrain","Hypothalamus","Hindbrain",
    "Amygdala","Midbrain","Hypothalamus","Medulla",
    "Hypothalamus","Hypothalamus","Thalamus","Pallidum/Septum",
    "Midbrain","Olfactory","Hindbrain","Pallidum/Septum",
    "Hindbrain","Midbrain"
  ),
  # FDR (q<0.05) status, read from the p_fdr column of
  # Drug_Statistical_Results_Primary.xlsx ("Veh vs Mor" sheet).
  # 21 of these top-30 |d| regions survive FDR; 9 are uncorrected p<0.05 only
  # (AVPV, ATN, PVH, RT, AT, TT, VeCB, GPe, PRM).
  fdr_sig = c(
    TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  FALSE, TRUE,  TRUE,
    FALSE, TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,
    FALSE, TRUE,  FALSE, TRUE,  FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
  )
) %>%
  mutate(
    label      = factor(label, levels = rev(label)),
    structure  = factor(unname(div_map[acronym]), levels = names(structure_colors)),
    point_fill = dplyr::case_when(
      fdr_sig %in% TRUE  ~ unname(structure_colors[as.character(structure)]),
      fdr_sig %in% FALSE ~ "white",
      TRUE               ~ "grey85"   # no FDR status recorded (NA)
    )
  )

# Guard: warn if any region has no FDR status recorded
.fig2a_missing_fdr <- fig2a_data$acronym[is.na(fig2a_data$fdr_sig)]
if (length(.fig2a_missing_fdr))
  message("Panel A - no FDR status recorded (rendered grey) for: ",
          paste(.fig2a_missing_fdr, collapse = ", "))

p_fig2a <- ggplot(fig2a_data, aes(x = abs_d, y = label, color = structure)) +
  
  geom_segment(aes(x = 0, xend = abs_d, yend = label),
               linewidth = 0.55, lineend = "round") +
  
  # Point fill encodes FDR status: filled (structure color) = q<0.05,
  # open (white) = uncorrected p<0.05 only, grey = no FDR status recorded.
  # Ring color (inherited from top aes) always shows structure.
  geom_point(aes(fill = point_fill,
                 shape = ifelse(fdr_sig %in% TRUE, "FDR q<0.05", "Uncorrected p<0.05")),
             size = 2.4, stroke = 0.6,
             show.legend = c(colour = FALSE, fill = FALSE, shape = TRUE)) +
  
  scale_color_manual(values = structure_colors, name = "Structure",
                     guide = guide_legend(order = 1)) +
  scale_fill_identity() +
  scale_shape_manual(
    name   = NULL,
    breaks = c("FDR q<0.05", "Uncorrected p<0.05"),
    values = c("FDR q<0.05" = 21, "Uncorrected p<0.05" = 21),
    guide  = guide_legend(order = 2,
                          override.aes = list(fill = c("grey35", "white"), colour = "grey35"))) +
  
  scale_x_continuous(
    name   = "Cohen's |d|",
    limits = c(0, 2.25),
    breaks = c(0, 0.5, 1.0, 1.5, 2.0),
    expand = c(0, 0)
  ) +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x        = element_text(size = 8, margin = margin(t = 4)),
    axis.title.y        = element_blank(),
    axis.text.x         = element_text(size = 7, color = "black"),
    axis.text.y         = element_text(size = 7, color = "black", hjust = 1),
    axis.line.x         = element_line(linewidth = 0.4, color = "#333333"),
    axis.line.y         = element_blank(),
    axis.ticks.y        = element_blank(),
    axis.ticks.x        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.grid.major.x  = element_line(linewidth = 0.2, color = "#eeeeee"),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.title        = element_text(size = 6, face = "bold"),
    legend.text         = element_text(size = 6),
    legend.key.size     = unit(8, "pt"),
    legend.spacing.y    = unit(2, "pt"),
    legend.position     = "right",
    legend.box.spacing  = unit(2, "pt"),   # pull legend in close to the plot
    legend.background   = element_rect(fill = "white", color = NA),
    legend.margin       = margin(2, 2, 2, 2),
    plot.margin         = margin(8, 2, 6, 4),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.caption        = element_text(size = 6, color = "#555555", hjust = 0)
  ) +
  
  labs(tag = NULL)

save_panel(p_fig2a, file.path(output_dir, "Figure2_PanelA.tiff"),
           width = 114, height = 110)   # composite slot (Fig 2 left column)
cat("Panel A saved (lollipop).\n")


# =============================================================================
# PANEL E - Faceted exemplar bars: MBO, ASO, TRS, MEA; all 4 conditions; combined sex
# =============================================================================

# -- Load CCFv3 volumes --------------------------------------------------------
vol_raw <- read_excel(vol_file, skip = 1)
vol_map <- setNames(
  as.numeric(vol_raw[["Mean Volume (m)"]]),
  trimws(as.character(vol_raw[["abbreviation"]]))
)

# -- Load and normalize per-animal data ---------------------------------------
load_and_normalize <- function(region_acronym, sheet_name, sex_label) {
  raw <- read_excel(raw_file, sheet = sheet_name)
  
  acr_col    <- grep("acronym", names(raw), ignore.case = TRUE)[1]
  region_row <- which(trimws(as.character(raw[[acr_col]])) == region_acronym)
  if (length(region_row) == 0)
    stop(paste("Region not found:", region_acronym, "in", sheet_name))
  
  region_vol <- vol_map[region_acronym]
  if (is.na(region_vol) || region_vol == 0)
    stop(paste("Volume not found:", region_acronym))
  
  all_cols  <- names(raw)
  data_cols <- grep("^(Vehicle|Morphine|Chronic|Ro)", all_cols, value = TRUE)
  
  cond_map <- list(
    "Vehicle"            = grep("^Vehicle",  data_cols, value = TRUE),
    "Acute Morphine"     = grep("^Morphine", data_cols, value = TRUE),
    "Morphine-Dependent" = grep("^Chronic",  data_cols, value = TRUE),
    "Ro 64-6198"         = grep("^Ro",       data_cols, value = TRUE)
  )
  
  rows_list <- list()
  for (cond in names(cond_map)) {
    cols <- cond_map[[cond]]
    if (length(cols) == 0) next
    vals <- suppressWarnings(as.numeric(raw[region_row, cols]))
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals) > 0) {
      rows_list[[cond]] <- tibble(
        condition = cond,
        sex       = sex_label,
        value     = log10(vals / region_vol)
      )
    }
  }
  bind_rows(rows_list)
}

# -- Significance lookup from statistical result files ------------------------
# Reads sig_fdr (Benjamini-Hochberg corrected) directly from the ANOVA output so
# brackets always reflect the current statistics and are consistent with the FDR
# encoding used in the lollipop (Fig 2) and diverging chart (Fig 3A). Rerunning
# the analysis flows straight through here.
stat_primary_file <- file.path(OUTPUT_DIR, "Drug_Statistical_Results_Primary.xlsx")
stat_bysex_file   <- file.path(OUTPUT_DIR, "Drug_Statistical_Results_BySex.xlsx")

# Map comparison name -> x-axis positions (Veh=1, Mor=2, MorDep=3, Ro=4)
comparison_positions <- list(
  "Veh vs Mor"    = c(1, 2),
  "Veh vs MorDep" = c(1, 3),
  "Veh vs Ro"     = c(1, 4),
  "Mor vs MorDep" = c(2, 3),
  "Mor vs Ro"     = c(2, 4)
)

# Look up sig_fdr for one region in one sheet; returns "ns" if absent
get_sig <- function(stat_file, sheet, region_acronym) {
  df <- suppressMessages(read_excel(stat_file, sheet = sheet))
  acr_col <- grep("acronym", names(df), ignore.case = TRUE)[1]
  row <- df[trimws(as.character(df[[acr_col]])) == region_acronym, ]
  if (nrow(row) == 0) return("ns")
  sig <- as.character(row[["sig_fdr"]][1])
  if (is.na(sig)) return("ns")
  sig
}

# Build bracket list (combined sex - Primary file). Comparisons listed narrowest
# first so they stack without overlapping. Only FDR q<0.05 and below (*, **, ***)
# produce a bracket; FDR trends q<0.10 (dagger) and ns are not drawn.
build_brackets_combined <- function(region_acronym, comparisons) {
  brackets <- list()
  for (comp in comparisons) {
    sig <- get_sig(stat_primary_file, comp, region_acronym)
    if (sig %in% c("*", "**", "***")) {
      pos <- comparison_positions[[comp]]
      brackets[[length(brackets) + 1]] <- list(x1 = pos[1], x2 = pos[2], sig = sig)
    }
  }
  brackets
}

# Build bracket list within a sex (BySex file). Same p<0.05 threshold.
build_brackets_sex <- function(region_acronym, sex_label, comparisons) {
  sex_suffix <- if (sex_label == "Males") "Males" else "Fem"
  brackets <- list()
  for (comp in comparisons) {
    sheet <- paste0(comp, " - ", sex_suffix)
    sig <- get_sig(stat_bysex_file, sheet, region_acronym)
    if (sig %in% c("*", "**", "***")) {
      pos <- comparison_positions[[comp]]
      brackets[[length(brackets) + 1]] <- list(x1 = pos[1], x2 = pos[2], sig = sig)
    }
  }
  brackets
}

# Region full names for titles
region_labels <- c(
  MBO = "Mammillary body",
  ASO = "Accessory supraoptic group",
  TRS = "Triangular nucleus of septum",
  MEA = "Medial amygdalar nucleus"
)

# Significance annotations - read from ANOVA output at runtime (vs Vehicle)
# Combined sex: Veh vs Mor (1-2), Veh vs MorDep (1-3), Veh vs Ro (1-4)
fig2d_comparisons <- c("Veh vs Mor", "Veh vs MorDep", "Veh vs Ro")

# -- Add brackets via annotate() -----------------------------------------------
add_brackets <- function(p, brackets, y_data_max, y_min) {
  step   <- (y_data_max - y_min) * 0.14
  step   <- max(step, 0.13)
  tick   <- 0.05
  sig_only <- Filter(function(b) b$sig != "ns", brackets)
  for (i in seq_along(sig_only)) {
    b <- sig_only[[i]]
    y <- y_data_max + 0.10 + step * (i - 1)
    p <- p +
      annotate("segment", x=b$x1, xend=b$x2, y=y, yend=y,
               color="black", linewidth=0.32) +
      annotate("segment", x=b$x1, xend=b$x1, y=y, yend=y-tick,
               color="black", linewidth=0.32) +
      annotate("segment", x=b$x2, xend=b$x2, y=y, yend=y-tick,
               color="black", linewidth=0.32) +
      annotate("text", x=(b$x1+b$x2)/2, y=y+0.05, label=b$sig,
               size=2.4, family=FONT_FAMILY, vjust=0)
  }
  p
}

# -- Faceted exemplar bars (one panel: MBO, ASO, TRS, MEA) ---------------------
# Single faceted panel with a shared x-axis and per-region (free) y-scales,
# replacing the four standalone plots. This removes the repeated y-axis and
# x-axis ink and reads as one unit. Significance brackets are added per facet
# from the same FDR lookup used everywhere else.
cat("Loading exemplar-bar data...\n")
bar_regions <- c("MBO","ASO","TRS","MEA")

bar_df <- purrr::map_dfr(bar_regions, function(acr) {
  m <- load_and_normalize(acr, "Clustered Males",   "Males")
  f <- load_and_normalize(acr, "Clustered Females", "Females")
  bind_rows(m, f) %>% dplyr::select(-sex) %>% dplyr::filter(value > 0) %>%
    dplyr::mutate(region = acr)
}) %>%
  dplyr::mutate(condition = factor(condition, levels = condition_order),
                region    = factor(region,    levels = bar_regions))

# Per-region y-limits and bracket positions (facets use free y-scales)
brk_list <- list(); lim_list <- list()
for (acr in bar_regions) {
  d    <- dplyr::filter(bar_df, region == acr)
  ymax <- max(d$value, na.rm = TRUE); ymin <- min(d$value, na.rm = TRUE)
  brks <- Filter(function(b) b$sig %in% c("*","**","***"),
                 build_brackets_combined(acr, fig2d_comparisons))
  step <- max((ymax - ymin) * 0.14, 0.13)
  if (length(brks)) {
    for (i in seq_along(brks)) {
      b <- brks[[i]]
      brk_list[[length(brk_list) + 1]] <- data.frame(
        region = acr, x1 = b$x1, x2 = b$x2, xm = (b$x1 + b$x2) / 2,
        y = ymax + 0.10 + step * (i - 1), sig = b$sig,
        stringsAsFactors = FALSE)
    }
    ytop <- ymax + 0.10 + step * (length(brks) - 1) + step * 0.9
  } else {
    ytop <- ymax + 0.15
  }
  lim_list[[length(lim_list) + 1]] <- data.frame(
    region = acr, ylo = ymin - (ymax - ymin) * 0.05, yhi = ytop)
}
brk_df <- if (length(brk_list)) do.call(rbind, brk_list) else
  data.frame(region = character(), x1 = double(), x2 = double(),
             xm = double(), y = double(), sig = character())
lim_df <- do.call(rbind, lim_list)
brk_df$region <- factor(brk_df$region, levels = bar_regions)
lim_df$region <- factor(lim_df$region, levels = bar_regions)

# Invisible points pin each facet's y-range; tick marks for the brackets
lim_long <- tidyr::pivot_longer(lim_df, c(ylo, yhi), values_to = "value") %>%
  dplyr::mutate(condition = factor("Vehicle", levels = condition_order))
brk_tick <- if (nrow(brk_df)) {
  tidyr::pivot_longer(brk_df, c(x1, x2), values_to = "x") %>%
    dplyr::mutate(yend = y - 0.05)
} else brk_df

p_bars <- ggplot(bar_df, aes(x = condition, y = value)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.45, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.45) +
  
  geom_blank(data = lim_long, aes(x = condition, y = value)) +
  
  geom_segment(data = brk_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = brk_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = brk_df, aes(x = xm, y = y + 0.04, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = 2.4, vjust = 0) +
  
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = c("Vehicle"="Veh","Acute Morphine"="Mor",
                              "Morphine-Dependent"="MorDep","Ro 64-6198"="Ro")) +
  scale_y_continuous(name = expression(log[10](cFos^"+" ~ cells/mm^3)),
                     breaks = pretty_breaks(n = 4),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
  facet_wrap(~ region, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    strip.text        = element_text(size = 8, face = "bold", margin = margin(b = 2)),
    strip.background  = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 8, margin = margin(r = 2)),
    axis.text.x       = element_text(size = 7, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    panel.spacing     = unit(5, "pt"),
    legend.position   = "none",
    plot.margin       = margin(4, 6, 2, 2),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_bars, file.path(output_dir, "Figure2_PanelE_bars.tiff"),
           width = 169, height = 46)   # composite slot (Fig 2 bottom row, full width)
cat("Panel E (faceted exemplar bars) saved.\n")


# =============================================================================
# PANEL B - Significant-region count summary (combined sex)
# Number of regions significant vs. Vehicle for each drug condition, at two
# thresholds: uncorrected p<0.05 (faded) and FDR q<0.05 (solid, a subset of the
# uncorrected count). Values from Drug_Statistical_Results_Primary.xlsx
# (Vehicle < condition elevations): Acute 74/26, Dependent 50/18,
# Ro 3/0. The faded->solid collapse at Ro shows the restricted NOP-agonist
# profile (3 regions, none surviving FDR). Full counts in Table S2.
# =============================================================================

count_levels <- c("Acute Morphine", "Morphine-Dependent", "Ro 64-6198")
count_data <- tibble(
  condition = factor(rep(count_levels, 2), levels = count_levels),
  tier      = factor(rep(c("Uncorrected p<0.05", "FDR q<0.05"), each = 3),
                     levels = c("Uncorrected p<0.05", "FDR q<0.05")),
  count     = c(74, 50, 3,   26, 18, 0)
) %>% arrange(tier)   # uncorrected drawn first (behind); FDR overlaid on top

p_counts <- ggplot(count_data, aes(condition, count, fill = condition, alpha = tier)) +
  
  geom_col(position = "identity", width = 0.68) +
  
  # Count labels: uncorrected total above each bar, FDR count atop its solid part
  geom_text(data = subset(count_data, tier == "Uncorrected p<0.05"),
            aes(label = count), vjust = -0.4, size = 2.8, fontface = "bold",
            color = "#555555", family = FONT_FAMILY, show.legend = FALSE) +
  geom_text(data = subset(count_data, tier == "FDR q<0.05" & count > 0),
            aes(label = count), vjust = -0.4, size = 2.8, fontface = "bold",
            color = "black", family = FONT_FAMILY, show.legend = FALSE) +
  
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_alpha_manual(values = c("Uncorrected p<0.05" = 0.32, "FDR q<0.05" = 1.0),
                     name = NULL) +
  
  scale_y_continuous(name = "Significant regions (vs. Vehicle)",
                     limits = c(0, 82), breaks = seq(0, 80, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = c("Acute Morphine"     = "Mor",
                              "Morphine-Dependent" = "MorDep",
                              "Ro 64-6198"         = "Ro")) +
  
  guides(alpha = guide_legend(override.aes = list(fill = "#777777"))) +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x         = element_blank(),
    axis.title.y         = element_text(size = 8, margin = margin(r = 3)),
    axis.text.x          = element_text(size = 7, color = "black",
                                        angle = 30, hjust = 1, vjust = 1),
    axis.text.y          = element_text(size = 7, color = "black"),
    axis.line            = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks           = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length    = unit(2, "pt"),
    legend.position      = "top",              # out of the plot -> no overlap with bars
    legend.direction     = "horizontal",    # single row, right-aligned
    legend.justification = "right",
    legend.text          = element_text(size = 6),
    legend.key.size      = unit(7, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.margin        = margin(0, 0, 1, 0),
    legend.background    = element_rect(fill = "white", color = NA),
    plot.margin          = margin(4, 4, 4, 4),
    plot.background      = element_rect(fill = "white", color = NA)
  )

save_panel(p_counts, file.path(output_dir, "Figure2_PanelB_CountSummary.tiff"),
           width = 50, height = 47)   # composite slot (Fig 2 right sidebar, top)
cat("Panel B (count summary) saved.\n")


# =============================================================================
# PANEL C - Whole-brain total cFos+ cells per animal, by condition (combined sex)
# Per-animal total = sum of per-region counts across all 198 analysis regions,
# pooled across sex. Supports the non-significant whole-brain total: a one-way
# ANOVA omnibus p is computed at runtime and annotated as a bare p-value (no
# "n.s." marker). Raw counts (zeros retained - they are real here, unlike in the
# log-density panels, where log10(0) is dropped).
# NOTE: assumes the 198 rows are the non-overlapping analysis region set (no
# parent/child double-counting); the absolute total scales with that set but the
# across-condition comparison is unaffected.
# =============================================================================

read_totals <- function(sheet, sex_label) {
  raw       <- read_excel(raw_file, sheet = sheet)
  data_cols <- grep("^(Vehicle|Morphine|Chronic|Ro)", names(raw), value = TRUE)
  cond_map  <- list(
    "Vehicle"            = grep("^Vehicle",  data_cols, value = TRUE),
    "Acute Morphine"     = grep("^Morphine", data_cols, value = TRUE),
    "Morphine-Dependent" = grep("^Chronic",  data_cols, value = TRUE),
    "Ro 64-6198"         = grep("^Ro",       data_cols, value = TRUE)
  )
  rows_list <- list()
  for (cond in names(cond_map)) {
    cols <- cond_map[[cond]]
    if (length(cols) == 0) next
    totals <- vapply(cols, function(cc)
      sum(suppressWarnings(as.numeric(raw[[cc]])), na.rm = TRUE), numeric(1))
    rows_list[[cond]] <- tibble(condition = cond, sex = sex_label,
                                animal = cols, total = as.numeric(totals))
  }
  bind_rows(rows_list)
}

total_data <- bind_rows(
  read_totals("Clustered Males",   "Males"),
  read_totals("Clustered Females", "Females")
) %>%
  mutate(condition = factor(condition, levels = condition_order))

# Omnibus one-way ANOVA on whole-brain totals (combined sex)
.aov_e <- aov(total ~ condition, data = total_data)
.p_e   <- summary(.aov_e)[[1]][["Pr(>F)"]][1]
.p_lab <- sprintf("ANOVA p = %.2f", .p_e)
cat(sprintf("Panel C - whole-brain total ANOVA: p = %.4f\n", .p_e))

p_e <- ggplot(total_data, aes(x = condition, y = total)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.5) +
  
  scale_color_manual(values = condition_colors, guide = "none") +
  
  scale_y_continuous(
    name   = expression("Total cFos"^"+" ~ "cells (whole brain)"),
    labels = scales::comma,
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  
  scale_x_discrete(
    labels = c("Vehicle"            = "Veh",
               "Acute Morphine"     = "Mor",
               "Morphine-Dependent" = "MorDep",
               "Ro 64-6198"         = "Ro")
  ) +
  
  annotate("text", x = 2.5, y = Inf, vjust = 1.6, label = .p_lab,
           size = 2.4, family = FONT_FAMILY, color = "#555555") +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 8, margin = margin(r = 3)),
    axis.text.x       = element_text(size = 7, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    legend.position   = "none",
    plot.margin       = margin(6, 8, 4, 4),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_e, file.path(output_dir, "Figure2_PanelC_TotalCells.tiff"),
           width = 50, height = 58)   # composite slot (Fig 2 right sidebar, bottom)
cat("Panel C (whole-brain total cFos+ cells) saved.\n")


cat("\nAll Figure 2 panels saved to:", output_dir, "\n")


# =============================================================================
# FIGURE 3 - redirect output to Figure3 subdirectory
# =============================================================================
output_dir <- file.path(OUTPUT_DIR, "figures", "Figure3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PANEL A - Diverging bar chart: males vs. females, Veh vs. Acute Morphine
# Males extend right (+), females extend left (-)
#
# INCLUSION RULE (principled, stated): a region is plotted if it survives FDR
# (q<0.05) in EITHER sex, OR is significant at uncorrected p<0.05 in FEMALES.
# Under acute morphine no female region survives FDR, so this resolves to the
# 20 male FDR-significant regions + the 6 additional female uncorrected-p<0.05
# regions (the female footprint) = 26 regions. Values are read from
# Drug_Statistical_Results_BySex.xlsx ("Veh vs Mor - Males" / "- Fem").
#
# Opacity: full = FDR q<0.05, mid = uncorrected p<0.05 only, faint = ns.
# Tip dot = FDR-significant. Male-uncorrected-only regions (PAG, PPN, RN, DG,
# CA, GPi, etc.) are intentionally excluded by the rule; they remain in Table S2.
# =============================================================================

panel_a_data <- tribble(
  ~acronym, ~label,                                  ~d_m,  ~d_f,  ~sig_m, ~sig_f, ~fdr_m, ~fdr_f, ~structure,
  # ---- FDR-significant in males (20) ----
  "NB",   "Nucleus of brachium of inf. colliculus",  3.320, 0.754,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "SNr",  "Substantia nigra, reticular part",        3.273, 0.505,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "PRE",  "Presubiculum",                            3.220, 0.083,  TRUE,  FALSE,  TRUE,  FALSE, "Hippocampus",
  "TRS",  "Triangular nucleus of septum",            3.036, 0.883,  TRUE,  FALSE,  TRUE,  FALSE, "Pallidum/Septum",
  "SAG",  "Nucleus sagulum",                         2.754, 0.146,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "SCm",  "Superior colliculus, motor",              2.634, 0.124,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "MRN",  "Midbrain reticular nucleus",              2.546, 0.123,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "CP",   "Caudoputamen",                            2.421, 0.352,  TRUE,  FALSE,  TRUE,  FALSE, "Striatum",
  "MBO",  "Mammillary body",                         2.287, 1.810,  TRUE,  TRUE,   TRUE,  FALSE, "Hypothalamus",
  "POST", "Postsubiculum",                           2.269, 0.151,  TRUE,  FALSE,  TRUE,  FALSE, "Hippocampus",
  "BAC",  "Bed nucleus of anterior commissure",      2.240, 0.406,  TRUE,  FALSE,  TRUE,  FALSE, "Pallidum/Septum",
  "SCs",  "Superior colliculus, sensory",            2.233, 0.191,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "GPe",  "Globus pallidus, external",               2.189, 0.456,  TRUE,  FALSE,  TRUE,  FALSE, "Pallidum/Septum",
  "AT",   "Anterior tegmental nucleus",              2.123, 0.572,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "PGRN", "Paragigantocellular reticular nucleus",   2.040, 0.794,  TRUE,  FALSE,  TRUE,  FALSE, "Medulla",
  "RAmb", "Midbrain raphe nuclei",                   1.943, 0.992,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "VTA",  "Ventral tegmental area",                  1.841, 0.108,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "RR",   "Retrorubral area",                        1.840, 0.442,  TRUE,  FALSE,  TRUE,  FALSE, "Midbrain",
  "MEPO", "Median preoptic nucleus",                 1.826, 0.798,  TRUE,  FALSE,  TRUE,  FALSE, "Hypothalamus",
  "ASO",  "Accessory supraoptic group",              1.797, 2.055,  TRUE,  FALSE,  TRUE,  FALSE, "Hypothalamus",
  # ---- Female footprint: uncorrected p<0.05 in females, not male-FDR (6) ----
  "SBPV", "Subparaventricular zone",                 1.784, 1.779,  TRUE,  TRUE,   FALSE, FALSE, "Hypothalamus",
  "AHN",  "Anterior hypothalamic nucleus",           1.760, 1.505,  TRUE,  TRUE,   FALSE, FALSE, "Hypothalamus",
  "NDB",  "Diagonal band nucleus",                   1.725, 1.742,  TRUE,  TRUE,   FALSE, FALSE, "Pallidum/Septum",
  "LPO",  "Lateral preoptic area",                   1.497, 1.499,  TRUE,  TRUE,   FALSE, FALSE, "Hypothalamus",
  "AVPV", "Anteroventral periventricular nucleus",   1.491, 2.276,  FALSE, TRUE,   FALSE, FALSE, "Hypothalamus",
  "LC",   "Locus coeruleus",                         0.060, 1.542,  FALSE, TRUE,   FALSE, FALSE, "Hindbrain"
) %>%
  # Sort by male |d| descending - females mirror on left
  arrange(desc(d_m)) %>%
  mutate(
    label     = factor(label, levels = rev(label)),
    structure = factor(unname(div_map[acronym]), levels = names(structure_colors)),
    # Opacity tiers: full = FDR q<0.05, mid = uncorrected p<0.05 only, faint = ns.
    alpha_m   = dplyr::case_when(fdr_m ~ 1.0, sig_m ~ 0.55, TRUE ~ 0.22),
    alpha_f   = dplyr::case_when(fdr_f ~ 1.0, sig_f ~ 0.55, TRUE ~ 0.22)
  )

# Build diverging chart - males positive, females negative
# Reshape to long format for easier plotting
df_long_a <- bind_rows(
  panel_a_data %>%
    transmute(label, acronym, structure,
              value = d_m, side = "Males",
              alpha = alpha_m, sig = sig_m, fdr = fdr_m),
  panel_a_data %>%
    transmute(label, acronym, structure,
              value = -d_f, side = "Females",
              alpha = alpha_f, sig = sig_f, fdr = fdr_f)
) %>%
  dplyr::mutate(sig_tier = factor(dplyr::case_when(
    fdr ~ "FDR q<0.05", sig ~ "Uncorrected p<0.05", TRUE ~ "ns"),
    levels = c("FDR q<0.05", "Uncorrected p<0.05", "ns")))

p_a <- ggplot(df_long_a,
              aes(x = value, y = label, fill = structure, alpha = sig_tier)) +
  
  # Bars
  geom_col(width = 0.75, color = NA) +
  
  # Zero line
  geom_vline(xintercept = 0, linewidth = 0.5, color = "#333333") +
  
  # Solid dot at bar tip marks FDR-significant regions (q<0.05)
  geom_point(data = df_long_a %>% filter(fdr %in% TRUE),
             aes(x = value + sign(value) * 0.04, color = structure),
             size = 1.1, show.legend = FALSE) +
  
  scale_fill_manual(values  = structure_colors, name = "Structure",
                    guide = guide_legend(order = 1)) +
  scale_color_manual(values = structure_colors, name = "Structure", guide = "none") +
  scale_alpha_manual(
    name   = "vs. Vehicle",
    values = c("FDR q<0.05" = 1.0, "Uncorrected p<0.05" = 0.55, "ns" = 0.22),
    labels = c("FDR q<0.05  \u25CF", "Uncorrected p<0.05", "ns"),
    guide  = guide_legend(order = 2, override.aes = list(fill = "grey30"))) +
  
  scale_x_continuous(
    name   = "Cohen's |d|",
    limits = c(-2.6, 3.5),
    breaks = c(-2, -1, 0, 1, 2, 3),
    labels = c("2", "1", "0", "1", "2", "3"),
    expand = c(0, 0)
  ) +
  
  # Headroom above the top bar so the sex labels are not clipped
  scale_y_discrete(expand = expansion(add = c(0.6, 2.0))) +
  
  # Female / Male axis labels (placed in the headroom)
  annotate("text", x = -1.6, y = length(levels(df_long_a$label)) + 1.3,
           label = "\u2190 Females", hjust = 0.5,
           size = 3.0, family = FONT_FAMILY, fontface = "bold",
           color = sex_colors[["Females"]]) +
  annotate("text", x = 2.0, y = length(levels(df_long_a$label)) + 1.3,
           label = "Males \u2192", hjust = 0.5,
           size = 3.0, family = FONT_FAMILY, fontface = "bold",
           color = sex_colors[["Males"]]) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x   = element_text(size = 8, margin = margin(t = 4)),
    axis.title.y   = element_blank(),
    axis.text.x    = element_text(size = 7, color = "black"),
    axis.text.y    = element_text(size = 7, color = "black", hjust = 1),
    axis.line.x    = element_line(linewidth = 0.4, color = "#333333"),
    axis.line.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    axis.ticks.x   = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.grid.major.x  = element_line(linewidth = 0.2, color = "#eeeeee"),
    panel.grid.major.y  = element_blank(),
    legend.title        = element_text(size = 6, face = "bold"),
    legend.text         = element_text(size = 6),
    legend.key.size     = unit(8, "pt"),
    legend.position     = "right",
    legend.background   = element_rect(fill = "white", color = NA),
    plot.margin         = margin(18, 4, 6, 4),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.caption        = element_text(size = 6, color = "#555555", hjust = 0)
  ) +
  
  labs(tag = NULL)

save_panel(p_a, file.path(output_dir, "Figure3_PanelA.tiff"),
           width = 106, height = 122)   # composite slot (Fig 3 2x2, top-left, diverging bars)
cat("Panel A saved.\n")


# =============================================================================
# PANEL B - Sex-stratified bar graphs: AVPV, SNr, LC
# Males and females shown side-by-side within each region
# All 4 conditions; individual animals overlaid
# =============================================================================

# Load CCFv3 volumes

# Significance annotations for Panel B - read from ANOVA output at runtime.
# Within each sex: Veh vs Mor (1-2), Veh vs MorDep (1-3), Veh vs Ro (1-4),
# plus Mor vs MorDep (2-3). Brackets are drawn only where sig_fdr != "ns".
fig3b_comparisons <- c("Veh vs Mor", "Veh vs MorDep", "Veh vs Ro", "Mor vs MorDep")

# Add brackets helper. Draws significance stars (*, **, ***) only; comparisons
# that are not significant get no bracket.
add_brackets_sex <- function(p, brackets, y_max, y_min) {
  draw <- Filter(function(b) b$sig != "ns", brackets)
  if (length(draw) == 0) return(p)
  step <- max((y_max - y_min) * 0.14, 0.13)
  tick <- 0.05
  for (i in seq_along(draw)) {
    b <- draw[[i]]
    y <- y_max + 0.10 + step * (i - 1)
    p <- p +
      annotate("segment", x=b$x1, xend=b$x2, y=y, yend=y,
               color="black", linewidth=0.32) +
      annotate("segment", x=b$x1, xend=b$x1, y=y, yend=y-tick,
               color="black", linewidth=0.32) +
      annotate("segment", x=b$x2, xend=b$x2, y=y, yend=y-tick,
               color="black", linewidth=0.32) +
      annotate("text", x=(b$x1+b$x2)/2, y=y+0.05, label=b$sig,
               size=2.4, family=FONT_FAMILY, vjust=0)
  }
  p
}

# -- Faceted sex-stratified bars (one panel: region rows x sex columns) ---------
# Replaces the six standalone plots with a single faceted panel: SNr / AVPV / LC
# down the rows, Males / Females across the columns. Y-scales are free per region
# (shared across the two sexes within a region). Brackets are added per cell from
# the same within-sex FDR lookup.
cat("Loading Panel B data (faceted)...\n")
bsex_regions <- c("SNr","AVPV","LC")
sex_levels   <- c("Males","Females")

panelB_df <- purrr::map_dfr(bsex_regions, function(acr) {
  bind_rows(
    load_and_normalize(acr, "Clustered Males",   "Males"),
    load_and_normalize(acr, "Clustered Females", "Females")
  ) %>% dplyr::filter(value > 0) %>% dplyr::mutate(region = acr)
}) %>%
  dplyr::mutate(condition = factor(condition, levels = condition_order),
                region    = factor(region, levels = bsex_regions),
                sex       = factor(sex,    levels = sex_levels))

# Per-region y-limits (shared across sexes); per-(region,sex) bracket positions
brkB_list <- list(); limB_list <- list()
for (acr in bsex_regions) {
  d_reg <- dplyr::filter(panelB_df, region == acr)
  ymax  <- max(d_reg$value, na.rm = TRUE); ymin <- min(d_reg$value, na.rm = TRUE)
  step  <- max((ymax - ymin) * 0.14, 0.13)
  n_sig_max <- 0
  for (sx in sex_levels) {
    brks <- Filter(function(b) b$sig != "ns",
                   build_brackets_sex(acr, sx, fig3b_comparisons))
    if (length(brks)) {
      for (i in seq_along(brks)) {
        b <- brks[[i]]
        brkB_list[[length(brkB_list) + 1]] <- data.frame(
          region = acr, sex = sx, x1 = b$x1, x2 = b$x2,
          xm = (b$x1 + b$x2) / 2, y = ymax + 0.10 + step * (i - 1),
          sig = b$sig, stringsAsFactors = FALSE)
      }
      n_sig_max <- max(n_sig_max, length(brks))
    }
  }
  ytop <- if (n_sig_max > 0) ymax + 0.10 + step * (n_sig_max - 1) + step * 0.85
  else ymax + 0.15
  limB_list[[length(limB_list) + 1]] <- data.frame(
    region = acr, ylo = ymin - (ymax - ymin) * 0.05, yhi = ytop)
}
brkB_df <- if (length(brkB_list)) do.call(rbind, brkB_list) else
  data.frame(region = character(), sex = character(), x1 = double(),
             x2 = double(), xm = double(), y = double(), sig = character())
limB_df <- do.call(rbind, limB_list)
brkB_df$region <- factor(brkB_df$region, levels = bsex_regions)
limB_df$region <- factor(limB_df$region, levels = bsex_regions)
if (nrow(brkB_df)) brkB_df$sex <- factor(brkB_df$sex, levels = sex_levels)

# Pin each region row's y-range in both sex columns; bracket tick marks
limB_long <- limB_df %>%
  tidyr::pivot_longer(c(ylo, yhi), values_to = "value") %>%
  tidyr::crossing(sex = factor(sex_levels, levels = sex_levels)) %>%
  dplyr::mutate(condition = factor("Vehicle", levels = condition_order))
brkB_tick <- if (nrow(brkB_df)) {
  tidyr::pivot_longer(brkB_df, c(x1, x2), values_to = "x") %>%
    dplyr::mutate(yend = y - 0.05)
} else brkB_df

# Per-facet x-axis baseline - facet_grid + theme_classic drops the axis line on
# interior facets, so draw it at each facet's lower limit in both sex columns.
baseB <- limB_df %>%
  tidyr::crossing(sex = factor(sex_levels, levels = sex_levels))

p_bsex <- ggplot(panelB_df, aes(x = condition, y = value)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.4, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.4) +
  
  geom_blank(data = limB_long, aes(x = condition, y = value)) +
  
  geom_segment(data = baseB, aes(x = 0.5, xend = 4.5, y = ylo, yend = ylo),
               inherit.aes = FALSE, linewidth = 0.4, color = "#333333") +
  
  geom_segment(data = brkB_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_segment(data = brkB_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_text(data = brkB_df, aes(x = xm, y = y + 0.04, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = 2.4, vjust = 0) +
  
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = c("Vehicle"="Veh","Acute Morphine"="Mor",
                              "Morphine-Dependent"="MorDep","Ro 64-6198"="Ro")) +
  scale_y_continuous(name = expression(log[10](cFos^"+" ~ cells/mm^3)),
                     breaks = pretty_breaks(n = 4),
                     expand = expansion(mult = c(0, 0.04))) +
  
  facet_grid(region ~ sex, scales = "free_y") +
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    strip.text.x      = element_text(size = 8, face = "bold", margin = margin(b = 2)),
    strip.text.y      = element_text(size = 8, face = "bold", margin = margin(l = 2)),
    strip.background  = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 8, margin = margin(r = 2)),
    axis.text.x       = element_text(size = 7, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7, color = "black"),
    axis.line.y       = element_line(linewidth = 0.4, color = "#333333"),
    axis.line.x       = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    panel.spacing.x   = unit(4, "pt"),
    panel.spacing.y   = unit(5, "pt"),
    legend.position   = "none",
    plot.margin       = margin(4, 6, 2, 2),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_bsex, file.path(output_dir, "Figure3_PanelB_SexBars.tiff"),
           width = 58, height = 122)   # composite slot (Fig 3 2x2, top-right sidebar, sex bars (filled))
cat("Panel B (faceted sex bars) saved.\n")


# =============================================================================
# PANEL C - Acute vs. dependent effect-size comparison (combined sex)
# x = Cohen's |d|, Vehicle vs. acute morphine; y = Cohen's |d|, Vehicle vs.
# morphine-dependent. Each point is a region FDR-significant (q<0.05) in at
# least one of the two conditions (30 regions). Dashed line is y = x (equal
# effect size in both states).
#
# DESCRIPTIVE grouping by FDR status vs. VEHICLE only (NOT a test of change):
#   "in both"          = FDR vs vehicle in acute AND dependence (n=14)
#   "dependence only"  = FDR vs vehicle in dependence, not acute (n=4)
#   "acute only"       = FDR vs vehicle in acute, not dependence (n=12)
# IMPORTANT: the direct acute-vs-dependent contrast (Mor vs MorDep) yields
# 0/198 FDR-significant regions, so no point here differs significantly between
# states. Whole-brain effect sizes are moderately correlated (Pearson r=0.44,
# Spearman rho=0.33, all 198 regions). Panel is descriptive, not inferential.
# All effects are activations (Vehicle < condition). Values from p_fdr/cohens_d
# columns of Drug_Statistical_Results_Primary.xlsx. Table S2.
# Uses ggrepel for labels if installed; falls back to geom_text otherwise.
# =============================================================================

transition_data <- tribble(
  ~acronym, ~label,                      ~d_acute, ~d_dep, ~category,      ~lab,
  "MBO", "Mammillary body", 2.070, 2.062, "Maintained", TRUE,
  "PVpo", "PVN, preoptic", 1.775, 1.147, "Attenuated", FALSE,
  "NOD", "Nodulus", 1.742, 0.276, "Attenuated", TRUE,
  "ASO", "Accessory supraoptic", 1.945, 1.734, "Maintained", TRUE,
  "LPO", "Lateral preoptic", 1.545, 1.731, "Maintained", FALSE,
  "TRS", "Triangular n. septum", 1.747, 1.704, "Maintained", FALSE,
  "NDB", "Diagonal band", 1.702, 1.256, "Attenuated", FALSE,
  "MEA", "Medial amygdala", 1.431, 1.606, "Maintained", FALSE,
  "DMH", "Dorsomedial hypothalamus", 1.244, 1.592, "Emergent", TRUE,
  "MPN", "Medial preoptic n.", 1.418, 1.582, "Maintained", FALSE,
  "SI", "Substantia innominata", 1.357, 1.525, "Maintained", FALSE,
  "PVa", "PVN, anterior", 1.523, 1.105, "Attenuated", FALSE,
  "BMA", "Basomedial amygdala", 1.250, 1.515, "Maintained", FALSE,
  "HATA", "Hippocampo-amyg. transition", 1.118, 1.487, "Emergent", TRUE,
  "PD", "Posterodorsal preoptic", 0.868, 1.487, "Emergent", TRUE,
  "SBPV", "Subparaventricular zone", 1.706, 1.479, "Maintained", FALSE,
  "RAmb", "Midbrain raphe", 1.464, 1.008, "Attenuated", FALSE,
  "PVi", "PVN, intermediate", 1.459, 1.023, "Attenuated", FALSE,
  "PMd", "Dorsal premammillary", 1.150, 1.451, "Emergent", TRUE,
  "AHN", "Anterior hypothalamic n.", 1.530, 1.450, "Maintained", FALSE,
  "MEPO", "Median preoptic n.", 1.275, 1.443, "Maintained", FALSE,
  "NLL", "N. lateral lemniscus", 1.439, 0.357, "Attenuated", TRUE,
  "LHA", "Lateral hypothalamic", 1.504, 1.421, "Maintained", FALSE,
  "SNr", "Substantia nigra retic.", 1.419, 0.821, "Attenuated", TRUE,
  "PGRN", "Paragigantocell. retic.", 1.398, 0.612, "Attenuated", FALSE,
  "MPO", "Medial preoptic area", 1.167, 1.344, "Maintained", FALSE,
  "PeF", "Perifornical n.", 1.380, 1.296, "Maintained", FALSE,
  "PN", "Paranigral n.", 1.281, 0.983, "Attenuated", FALSE,
  "PIR", "Piriform area", 1.254, 1.285, "Attenuated", FALSE,
  "AAA", "Anterior amygdalar a.", 1.126, 1.045, "Attenuated", FALSE,
) %>%
  # Internal codes Maintained/Emergent/Attenuated are remapped to neutral,
  # descriptive labels (FDR-vs-vehicle status only; no change is implied).
  mutate(category = factor(category,
                           levels = c("Maintained", "Emergent", "Attenuated"),
                           labels = c("in both", "dependence only", "acute only")))

transition_colors <- c(
  "in both"         = "#2b8cbe",   # blue   - FDR vs vehicle in acute AND dependence
  "dependence only" = "#e6550d",   # orange - FDR vs vehicle in dependence only
  "acute only"      = "#969696"    # grey   - FDR vs vehicle in acute only
)

# Label layer: ggrepel if available, plain geom_text fallback so the script
# never breaks on a missing package.
.lab_df <- subset(transition_data, lab)
label_layer <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(
    data = .lab_df, aes(label = acronym),
    size = 2.3, family = FONT_FAMILY, fontface = "bold",
    min.segment.length = 0, segment.size = 0.25, segment.color = "#999999",
    box.padding = 0.35, point.padding = 0.2, max.overlaps = Inf,
    show.legend = FALSE
  )
} else {
  geom_text(data = .lab_df, aes(label = acronym),
            size = 2.3, family = FONT_FAMILY, fontface = "bold",
            vjust = -0.7, show.legend = FALSE)
}

p_dep <- ggplot(transition_data, aes(x = d_acute, y = d_dep, color = category)) +
  
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.4, color = "#bbbbbb") +
  
  geom_point(size = 2.0, alpha = 0.9) +
  label_layer +
  
  annotate("text", x = 0.05, y = 2.15, hjust = 0, vjust = 1,
           size = 2.5, family = FONT_FAMILY, color = "#555555",
           label = "Effect-size r = 0.44 (all 198 regions)") +
  annotate("text", x = 0.05, y = 1.98, hjust = 0, vjust = 1,
           size = 2.5, family = FONT_FAMILY, color = "#555555", fontface = "italic",
           label = "Direct acute-vs-dependent test: 0/198 FDR-sig") +
  
  scale_color_manual(values = transition_colors, name = "FDR-significant\nvs. vehicle:") +
  
  scale_x_continuous(name = expression("Acute Morphine  Cohen's |" * italic(d) * "|"),
                     limits = c(0, 2.2), breaks = seq(0, 2, 0.5), expand = c(0, 0)) +
  scale_y_continuous(name = expression("Morphine-Dependent  Cohen's |" * italic(d) * "|"),
                     limits = c(0, 2.2), breaks = seq(0, 2, 0.5), expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x      = element_text(size = 8, margin = margin(t = 4)),
    axis.title.y      = element_text(size = 8, margin = margin(r = 4)),
    axis.text         = element_text(size = 7, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    legend.position      = c(0.02, 0.52),   # inside, tucked against the y-axis (empty left region)
    legend.justification = c(0, 0.5),
    legend.title         = element_text(size = 6, face = "bold"),
    legend.text          = element_text(size = 6),
    legend.key.size   = unit(9, "pt"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin       = margin(6, 8, 6, 6),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_dep, file.path(output_dir, "Figure3_PanelC_AcuteVsDependence.tiff"),
           width = 106, height = 86)   # composite slot (Fig 3 2x2, bottom-left, scatter)
cat("Acute-vs-dependent effect-size comparison (Panel C) saved.\n")

# =============================================================================
# PANEL D - Acute vs. dependent exemplar bars (combined sex)
# Two regions embodying the off-diagonal categories from Panel C:
#   DMH = dependence-only (up in MorDep, not Mor)
#   NOD = acute-only      (up in Mor, not MorDep)
# Same faceted recipe / styling as the Figure 2 exemplar bars.
# =============================================================================
cat("Loading Panel D data (acute-vs-dependent exemplars)...\n")
pd_regions <- c("DMH","NOD")

pd_df <- purrr::map_dfr(pd_regions, function(acr) {
  m <- load_and_normalize(acr, "Clustered Males",   "Males")
  f <- load_and_normalize(acr, "Clustered Females", "Females")
  bind_rows(m, f) %>% dplyr::select(-sex) %>% dplyr::filter(value > 0) %>%
    dplyr::mutate(region = acr)
}) %>%
  dplyr::mutate(condition = factor(condition, levels = condition_order),
                region    = factor(region,    levels = pd_regions))

brkD_list <- list(); limD_list <- list()
for (acr in pd_regions) {
  d    <- dplyr::filter(pd_df, region == acr)
  ymax <- max(d$value, na.rm = TRUE); ymin <- min(d$value, na.rm = TRUE)
  brks <- Filter(function(b) b$sig %in% c("*","**","***"),
                 build_brackets_combined(acr, fig2d_comparisons))
  step <- max((ymax - ymin) * 0.14, 0.13)
  if (length(brks)) {
    for (i in seq_along(brks)) {
      b <- brks[[i]]
      brkD_list[[length(brkD_list) + 1]] <- data.frame(
        region = acr, x1 = b$x1, x2 = b$x2, xm = (b$x1 + b$x2) / 2,
        y = ymax + 0.10 + step * (i - 1), sig = b$sig,
        stringsAsFactors = FALSE)
    }
    ytop <- ymax + 0.10 + step * (length(brks) - 1) + step * 0.9
  } else {
    ytop <- ymax + 0.15
  }
  limD_list[[length(limD_list) + 1]] <- data.frame(
    region = acr, ylo = ymin - (ymax - ymin) * 0.05, yhi = ytop)
}
brkD_df <- if (length(brkD_list)) do.call(rbind, brkD_list) else
  data.frame(region = character(), x1 = double(), x2 = double(),
             xm = double(), y = double(), sig = character())
limD_df <- do.call(rbind, limD_list)
brkD_df$region <- factor(brkD_df$region, levels = pd_regions)
limD_df$region <- factor(limD_df$region, levels = pd_regions)

limD_long <- tidyr::pivot_longer(limD_df, c(ylo, yhi), values_to = "value") %>%
  dplyr::mutate(condition = factor("Vehicle", levels = condition_order))
brkD_tick <- if (nrow(brkD_df)) {
  tidyr::pivot_longer(brkD_df, c(x1, x2), values_to = "x") %>%
    dplyr::mutate(yend = y - 0.05)
} else brkD_df

p_dbars <- ggplot(pd_df, aes(x = condition, y = value)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.45, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.45) +
  
  geom_blank(data = limD_long, aes(x = condition, y = value)) +
  
  geom_segment(data = brkD_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = brkD_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = brkD_df, aes(x = xm, y = y + 0.04, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = 2.4, vjust = 0) +
  
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = c("Vehicle"="Veh","Acute Morphine"="Mor",
                              "Morphine-Dependent"="MorDep","Ro 64-6198"="Ro")) +
  scale_y_continuous(name = expression(log[10](cFos^"+" ~ cells/mm^3)),
                     breaks = pretty_breaks(n = 4),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
  facet_wrap(~ region, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    strip.text        = element_text(size = 8, face = "bold", margin = margin(b = 2)),
    strip.background  = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 8, margin = margin(r = 2)),
    axis.text.x       = element_text(size = 7, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    panel.spacing     = unit(5, "pt"),
    legend.position   = "none",
    plot.margin       = margin(4, 6, 2, 2),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_dbars, file.path(output_dir, "Figure3_PanelD_AcuteVsDep_bars.tiff"),
           width = 58, height = 86)   # composite slot (Fig 3 2x2, bottom-right, exemplars (filled))
cat("Panel D (acute-vs-dependent exemplar bars) saved.\n")

cat("\nAll Figure 3 panels saved to:", output_dir, "\n")