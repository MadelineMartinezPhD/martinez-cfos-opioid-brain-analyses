# regional_analysis.R
# =============================================================================
# Whole-brain cFos regional activation analysis.
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
# Repository:   https://github.com/MadelineMartinezPhD/martinez-cfos-opioid-brain-analyses
# Archived:     https://doi.org/10.5281/zenodo.21502549
# License:      MIT (see LICENSE)
# Tested with:  R 4.4.3 (see sessionInfo.txt)
#
# Usage (from the repository root):
#   Rscript scripts/regional_analysis.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# ANALYSIS PRIORITY:
#   PRIMARY   - One-way ANOVA on pooled males + females (condition only, no sex term)
#               Output: Drug_Statistical_Results_Primary.xlsx
#
#   SECONDARY - Sex-stratified one-way ANOVA (males and females separately)
#               Output: Drug_Statistical_Results_BySex.xlsx
#
#   SUPPLEMENTARY - Two-way ANOVA (condition x sex) and Negative Binomial
#               Output: Drug_Statistical_Results_TwoWay.xlsx
#                       Drug_Statistical_Results_NegBin.xlsx
#
# DATA STRUCTURE:
#   Drug_Composite.xlsx: "Clustered Males" and "Clustered Females" sheets
#   198 CCFv3 regions, 4 conditions, n=5-6 per group per sex
#
# NORMALIZATION:
#   ANOVA methods: log10(cells/mm3) using CCFv3 atlas volumes (Wang et al. 2020)
#   NegBin method: raw integer counts with log(volume) as offset
#
# EFFECT DIRECTION CONVENTION:
#   For every comparison (A vs B), all signed quantities are reported as
#   (B - A) = treatment - reference, so a treatment-induced INCREASE is POSITIVE.
#   This applies to estimate, t_statistic / z_statistic, cohens_d, and (via
#   exp) IRR, and is uniform across the primary, by-sex, two-way, and NegBin
#   tables. It matches the signed Cohen's d convention used in Figure 3C.
#   fold_change is likewise treatment/reference (mean_B / mean_A). Tukey p-values
#   are two-sided and unaffected by direction. The `direction` column reports the
#   raw-count ordering (A > B or A < B).
#
# PAIRWISE COMPARISONS (5 of 6 - Morphine-Dependent vs Ro excluded):
#   Vehicle vs Morphine, Vehicle vs Morphine-Dependent, Vehicle vs Ro,
#   Morphine vs Morphine-Dependent, Morphine vs Ro
#
# Required packages:
#   install.packages(c("readxl", "openxlsx", "dplyr", "MASS", "emmeans"))

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("readxl", "openxlsx", "dplyr", "MASS", "emmeans")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(readxl); library(openxlsx); library(dplyr); library(MASS); library(emmeans)

# ---- Config ------------------------------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/Drug_cFos_Analysis_Combined.R
# To use a different layout, override the two directories below, either by
# editing them here or by setting the environment variables before running.
DATA_DIR    <- Sys.getenv("CFOS_DATA_DIR",   "data")
OUTPUT_DIR  <- Sys.getenv("CFOS_OUTPUT_DIR", "results")
RAW_FILE    <- file.path(DATA_DIR, "Drug_Composite.xlsx")
VOLUME_FILE <- file.path(DATA_DIR, "ccfv3_volumes.xlsx")
FDR_ALPHA   <- 0.05

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ---------------------------------------------------------------
cat("Loading data...\n")
df_m <- as.data.frame(read_excel(RAW_FILE, sheet = "Clustered Males"))
df_f <- as.data.frame(read_excel(RAW_FILE, sheet = "Clustered Females"))

get_data_cols <- function(df) {
  cols <- grep("Vehicle|Morphine|Ro", names(df), value = TRUE)
  cols[!grepl("^Unnamed", cols)]
}
data_cols_m <- get_data_cols(df_m)
data_cols_f <- get_data_cols(df_f)

make_groups <- function(cols) {
  list(
    Vehicle            = grep("^Vehicle",  cols, value = TRUE),
    Morphine           = grep("^Morphine", cols, value = TRUE),
    `Morphine-Dependent` = grep("^Chronic",  cols, value = TRUE),
    Ro                 = grep("^Ro",       cols, value = TRUE)
  )
}
groups_m   <- make_groups(data_cols_m)
groups_f   <- make_groups(data_cols_f)
cond_names <- names(groups_m)

cat(sprintf("  Males:   %d regions, %d animals\n", nrow(df_m), length(data_cols_m)))
cat(sprintf("  Females: %d regions, %d animals\n", nrow(df_f), length(data_cols_f)))

# ---- Column-integrity guard --------------------------------------------------
# A single header typo (e.g. "Morphne" instead of "Morphine") makes a real
# animal invisible to get_data_cols()/make_groups() and it is dropped silently,
# corrupting every statistic for that condition. Any header that looks like an
# animal column (contains "Male"/"Female" followed by a number) MUST be both
# captured and assigned to exactly one condition group, otherwise the run aborts with an error.
check_columns <- function(df, data_cols, groups, sex_label) {
  looks_data <- grep("(Male|Female)[[:space:]]*[0-9]", names(df), value = TRUE)
  grouped    <- unlist(groups, use.names = FALSE)
  dropped    <- setdiff(looks_data, data_cols)   # look like data but not captured
  ungrouped  <- setdiff(data_cols,  grouped)     # captured but no condition group
  if (length(dropped) || length(ungrouped)) {
    stop(sprintf(paste0(
      "[%s] column-integrity check FAILED.\n",
      "  Look like data but not captured: %s\n",
      "  Captured but not grouped:        %s\n",
      "  Fix the header typo (e.g. 'Morphne' -> 'Morphine') and re-run."),
      sex_label,
      if (length(dropped))   paste(dropped,   collapse = ", ") else "none",
      if (length(ungrouped)) paste(ungrouped, collapse = ", ") else "none"))
  }
  cat(sprintf("  [%s] integrity OK: %d columns grouped (%s)\n",
              sex_label, length(data_cols),
              paste(sprintf("%s=%d", names(groups), lengths(groups)), collapse = ", ")))
}
check_columns(df_m, data_cols_m, groups_m, "Males")
check_columns(df_f, data_cols_f, groups_f, "Females")

# ---- Volumes (Wang et al. 2020 CCFv3 Table S4) --------------------------------
cat("Loading volumes...\n")
vol_df          <- as.data.frame(read_excel(VOLUME_FILE, skip = 1))
vol_map         <- setNames(vol_df[["Mean Volume (m)"]], vol_df[["abbreviation"]])
df_m$volume_mm3 <- vol_map[df_m$acronym]
df_f$volume_mm3 <- vol_map[df_f$acronym]
cat(sprintf("  Matched %d/%d regions\n", sum(!is.na(df_m$volume_mm3)), nrow(df_m)))

# ---- Log10 density normalization ---------------------------------------------
log10_density <- function(df, data_cols) {
  out <- df[, c("name", "acronym")]
  for (col in data_cols) {
    vals <- df[[col]] / df$volume_mm3
    vals[vals == 0] <- NA
    out[[col]] <- log10(vals)
  }
  out
}
norm_m <- log10_density(df_m, data_cols_m)
norm_f <- log10_density(df_f, data_cols_f)

# Build pooled normalized data (males + females combined, same column names preserved)
norm_comb   <- df_m[, c("name", "acronym")]
groups_comb <- list()
df_comb     <- df_m[, c("name", "acronym", "volume_mm3")]
for (g in cond_names) {
  for (col in groups_m[[g]]) { norm_comb[[col]] <- norm_m[[col]]; df_comb[[col]] <- df_m[[col]] }
  for (col in groups_f[[g]]) { norm_comb[[col]] <- norm_f[[col]]; df_comb[[col]] <- df_f[[col]] }
  groups_comb[[g]] <- c(groups_m[[g]], groups_f[[g]])
}
cat("Normalization: log10(cells/mm3)\n")

# ---- Helpers -----------------------------------------------------------------
# 5 biologically relevant pairs - Morphine-Dependent vs Ro excluded
pairs <- list(
  c("Vehicle",  "Morphine"),
  c("Vehicle",  "Morphine-Dependent"),
  c("Vehicle",  "Ro"),
  c("Morphine", "Morphine-Dependent"),
  c("Morphine", "Ro")
)

# emmeans::contrast(method="pairwise") returns rows in combn(levels) order:
# (1,2),(1,3),(1,4),(2,3),(2,4),(3,4) for levels in the given order. We select
# the row for a pair by POSITION, not by matching the printed contrast label.
# Label matching with grepl is unsafe here because "Morphine" is a substring of
# "Morphine-Dependent", so a fallback search silently returns the wrong contrast.
pairwise_row <- function(ga, gb, levs) {
  ia <- match(ga, levs); ib <- match(gb, levs)
  if (is.na(ia) || is.na(ib)) return(NA_integer_)
  lo <- min(ia, ib); hi <- max(ia, ib)
  combos <- utils::combn(length(levs), 2)      # columns are (i, j), i < j
  which(combos[1, ] == lo & combos[2, ] == hi)
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ "\u2020",
    TRUE      ~ "ns"
  )
}

# Count regions below a significance threshold, read from the p-value column
# directly so the count does not depend on the significance-marker encoding.
n_below <- function(p, alpha) sum(p < alpha, na.rm = TRUE)

add_fdr <- function(res, p_col, fdr_col, star_col) {
  valid <- !is.na(res[[p_col]])
  res[[fdr_col]] <- NA
  res[[fdr_col]][valid] <- p.adjust(res[[p_col]][valid], method = "BH")
  res[[star_col]] <- sig_stars(res[[fdr_col]])
  res
}

write_excel <- function(results_list, filename, sig_col = "sig_tukey") {
  wb  <- createWorkbook()
  hdr <- createStyle(fontName="Arial", fontSize=10, textDecoration="bold",
                     border="Bottom", borderColour="#000000")
  dat <- createStyle(fontName="Arial", fontSize=10)
  sig <- createStyle(fontName="Arial", fontSize=10, fgFill="#FFF2CC")
  for (label in names(results_list)) {
    res   <- results_list[[label]]
    sname <- gsub("Morphine-Dependent", "MorDep", label)
    sname <- gsub("Morphine", "Mor", sname)
    sname <- gsub("Vehicle",  "Veh", sname)
    sname <- gsub("Females",  "Fem", sname)
    sname <- substr(sname, 1, 31)
    addWorksheet(wb, sname)
    writeData(wb, sname, res, startRow=1, startCol=1, rowNames=FALSE)
    addStyle(wb, sname, hdr, rows=1, cols=1:ncol(res), gridExpand=TRUE)
    setRowHeights(wb, sname, rows=1, heights=18)
    for (i in seq_len(nrow(res))) {
      sig_val <- res[[sig_col]][i]
      is_sig  <- isTRUE(!is.na(sig_val) && sig_val != "ns")
      addStyle(wb, sname, if (is_sig) sig else dat,
               rows=i+1, cols=1:ncol(res), gridExpand=TRUE)
    }
    setColWidths(wb, sname, cols=1:ncol(res), widths="auto")
    freezePane(wb, sname, firstActiveRow=2)
  }
  saveWorkbook(wb, filename, overwrite=TRUE)
  cat(sprintf("  Saved: %s\n", basename(filename)))
}

# Shared one-way ANOVA function - used by both primary and by-sex analyses
run_oneway <- function(norm_df, df_raw_sex, groups_sex, label_suffix = "") {
  out <- list()
  for (pair in pairs) {
    ga <- pair[1]; gb <- pair[2]
    label   <- if (nchar(label_suffix)) paste(ga, "vs", gb, "-", label_suffix)
    else paste(ga, "vs", gb)
    cols_a  <- groups_sex[[ga]]; cols_b <- groups_sex[[gb]]
    
    rows_list <- vector("list", nrow(norm_df))
    for (i in seq_len(nrow(norm_df))) {
      acr   <- norm_df$acronym[i]; rgn <- norm_df$name[i]
      va    <- as.numeric(norm_df[i, cols_a, drop=TRUE]);    va    <- va[!is.na(va)]
      vb    <- as.numeric(norm_df[i, cols_b, drop=TRUE]);    vb    <- vb[!is.na(vb)]
      raw_a <- as.numeric(df_raw_sex[i, cols_a, drop=TRUE]); raw_a <- raw_a[!is.na(raw_a)]
      raw_b <- as.numeric(df_raw_sex[i, cols_b, drop=TRUE]); raw_b <- raw_b[!is.na(raw_b)]
      mean_a_raw <- if (length(raw_a)) mean(raw_a) else NA
      mean_b_raw <- if (length(raw_b)) mean(raw_b) else NA
      
      all_v <- c(); all_l <- c()
      for (g in cond_names) {
        v <- as.numeric(norm_df[i, groups_sex[[g]], drop=TRUE]); v <- v[!is.na(v)]
        all_v <- c(all_v, v); all_l <- c(all_l, rep(g, length(v)))
      }
      
      empty <- data.frame(region_name=rgn, acronym=acr, comparison=paste(ga,"-",gb),
                          omnibus_F=NA, omnibus_p=NA, contrast_F=NA,
                          estimate=NA, std_error=NA, t_statistic=NA, p_tukey=NA,
                          fold_change=NA, cohens_d=NA,
                          group_a_mean=mean_a_raw, group_b_mean=mean_b_raw,
                          direction=NA, stringsAsFactors=FALSE)
      
      if (length(va) < 2 || length(vb) < 2) { rows_list[[i]] <- empty; next }
      
      tryCatch({
        # Explicit factor levels so Tukey row lookup is deterministic (never
        # dependent on R's default alphabetical ordering).
        fit   <- aov(all_v ~ factor(all_l, levels = cond_names))
        aov_s <- summary(fit)[[1]]
        ms_e  <- aov_s[["Mean Sq"]][2]
        tukey <- TukeyHSD(fit)[[1]]
        # p adj is two-sided (direction-independent); take it from whichever
        # row name matches this pair.
        key1  <- paste0(gb,"-",ga); key2 <- paste0(ga,"-",gb)
        p_tukey <- if (key1 %in% rownames(tukey)) tukey[key1, "p adj"] else
          if (key2 %in% rownames(tukey)) tukey[key2, "p adj"] else NA
        
        # All signed quantities: treatment - reference = (B - A), increase-positive.
        n_a <- length(va); n_b <- length(vb)
        estimate <- mean(vb) - mean(va)            # log10 mean difference (B - A)
        se   <- sqrt(ms_e * (1/n_a + 1/n_b))
        t_s  <- if (!is.na(estimate) && se > 0) estimate / se else NA
        fc   <- if (!is.na(mean_a_raw) && mean_a_raw > 0) mean_b_raw / mean_a_raw else NA
        sd_p <- sqrt(((n_a-1)*var(va) + (n_b-1)*var(vb)) / (n_a+n_b-2))
        cd   <- if (sd_p > 0) (mean(vb) - mean(va)) / sd_p else NA   # Cohen's d, B - A
        dir  <- if (!is.na(mean_a_raw) && mean_a_raw > mean_b_raw)
          paste(ga,">",gb) else paste(ga,"<",gb)
        
        rows_list[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                     comparison=paste(ga,"-",gb),
                                     omnibus_F=aov_s[["F value"]][1], omnibus_p=aov_s[["Pr(>F)"]][1],
                                     contrast_F=if(!is.na(t_s)) t_s^2 else NA,
                                     estimate=estimate, std_error=se, t_statistic=t_s, p_tukey=p_tukey,
                                     fold_change=fc, cohens_d=cd,
                                     group_a_mean=mean_a_raw, group_b_mean=mean_b_raw,
                                     direction=dir, stringsAsFactors=FALSE)
      }, error=function(e) { rows_list[[i]] <<- empty })
    }
    
    res <- do.call(rbind, rows_list)
    valid_t <- !is.na(res$p_tukey);  res$p_fdr <- NA
    res$p_fdr[valid_t]    <- p.adjust(res$p_tukey[valid_t], method="BH")
    valid_a <- !is.na(res$omnibus_p); res$omnibus_fdr <- NA
    res$omnibus_fdr[valid_a] <- p.adjust(res$omnibus_p[valid_a], method="BH")
    res$sig_tukey       <- sig_stars(res$p_tukey)
    res$sig_fdr         <- sig_stars(res$p_fdr)
    res$sig_omnibus_fdr <- sig_stars(res$omnibus_fdr)
    
    out[[label]] <- res
    cat(sprintf("    %-40s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
                n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
                n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
  }
  out
}

# ==============================================================================
# PRIMARY: One-way ANOVA - males + females pooled, condition only
# ==============================================================================
cat("\n--- PRIMARY: One-way ANOVA (pooled, no sex term) ---\n")
primary_results <- run_oneway(norm_comb, df_comb, groups_comb, label_suffix = "")
write_excel(primary_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_Primary.xlsx"))

# ==============================================================================
# SECONDARY: Sex-stratified one-way ANOVA
# ==============================================================================
cat("\n--- SECONDARY: Sex-stratified one-way ANOVA ---\n")
bysex_results <- c(
  run_oneway(norm_m, df_m, groups_m, label_suffix = "Males"),
  run_oneway(norm_f, df_f, groups_f, label_suffix = "Females")
)
write_excel(bysex_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_BySex.xlsx"))

# ==============================================================================
# SUPPLEMENTARY 1: Two-way ANOVA - condition x sex
# ==============================================================================
cat("\n--- SUPPLEMENTARY: Two-way ANOVA (condition x sex) ---\n")

make_long_twoway <- function(i) {
  rows <- list()
  for (g in cond_names) {
    for (col in groups_m[[g]]) {
      v <- norm_m[i, col]
      if (!is.na(v)) rows[[length(rows)+1]] <-
          data.frame(value=v, condition=g, sex="Male",   stringsAsFactors=FALSE)
    }
    for (col in groups_f[[g]]) {
      v <- norm_f[i, col]
      if (!is.na(v)) rows[[length(rows)+1]] <-
          data.frame(value=v, condition=g, sex="Female", stringsAsFactors=FALSE)
    }
  }
  long <- do.call(rbind, rows)
  long$condition <- factor(long$condition, levels=cond_names)
  long$sex       <- factor(long$sex, levels=c("Male","Female"))
  long
}

twoway_results <- list()
sex_rows       <- vector("list", nrow(norm_m))
int_rows       <- vector("list", nrow(norm_m))

for (i in seq_len(nrow(norm_m))) {
  acr  <- norm_m$acronym[i]; rgn <- norm_m$name[i]
  long <- make_long_twoway(i)
  sex_rows[[i]] <- data.frame(region_name=rgn, acronym=acr, F_sex=NA, p_sex=NA, stringsAsFactors=FALSE)
  int_rows[[i]] <- data.frame(region_name=rgn, acronym=acr, F_interaction=NA, p_interaction=NA, stringsAsFactors=FALSE)
  if (nrow(long) < 6) next
  tryCatch({
    fit    <- aov(value ~ condition * sex, data=long)
    summ   <- summary(fit)[[1]]
    rnames <- trimws(rownames(summ))
    sex_idx <- which(rnames == "sex")
    int_idx <- which(rnames == "condition:sex")
    sex_rows[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                F_sex = if (length(sex_idx)) summ[["F value"]][sex_idx] else NA,
                                p_sex = if (length(sex_idx)) summ[["Pr(>F)"]][sex_idx]  else NA,
                                stringsAsFactors=FALSE)
    int_rows[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                F_interaction = if (length(int_idx)) summ[["F value"]][int_idx] else NA,
                                p_interaction = if (length(int_idx)) summ[["Pr(>F)"]][int_idx]  else NA,
                                stringsAsFactors=FALSE)
  }, error=function(e) {
    sex_rows[[i]] <<- data.frame(region_name=rgn, acronym=acr, F_sex=NA, p_sex=NA, stringsAsFactors=FALSE)
    int_rows[[i]] <<- data.frame(region_name=rgn, acronym=acr, F_interaction=NA, p_interaction=NA, stringsAsFactors=FALSE)
  })
}

sex_effect_df <- add_fdr(do.call(rbind, sex_rows), "p_sex",         "q_sex",         "sig_sex")
int_effect_df <- add_fdr(do.call(rbind, int_rows), "p_interaction", "q_interaction", "sig_interaction")
cat(sprintf("  Sex main effect:               %3d (%3d) regions\n",
            n_below(sex_effect_df$q_sex, 0.05), n_below(sex_effect_df$q_sex, 0.10)))
cat(sprintf("  Condition x Sex interaction:   %3d (%3d) regions\n",
            n_below(int_effect_df$q_interaction, 0.05), n_below(int_effect_df$q_interaction, 0.10)))

for (pair in pairs) {
  ga <- pair[1]; gb <- pair[2]; label <- paste(ga, "vs", gb)
  rows_list <- vector("list", nrow(norm_m))
  
  for (i in seq_len(nrow(norm_m))) {
    acr  <- norm_m$acronym[i]; rgn <- norm_m$name[i]
    long <- make_long_twoway(i)
    raw_a <- c(as.numeric(df_m[i, groups_m[[ga]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[ga]], drop=TRUE]))
    raw_b <- c(as.numeric(df_m[i, groups_m[[gb]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[gb]], drop=TRUE]))
    raw_a <- raw_a[!is.na(raw_a)]; raw_b <- raw_b[!is.na(raw_b)]
    mean_a <- if (length(raw_a)) mean(raw_a) else NA
    mean_b <- if (length(raw_b)) mean(raw_b) else NA
    
    empty <- data.frame(region_name=rgn, acronym=acr, comparison=paste(ga,"-",gb),
                        omnibus_F=NA, omnibus_p=NA, contrast_F=NA,
                        estimate=NA, std_error=NA, t_statistic=NA, p_tukey=NA,
                        fold_change=NA, cohens_d=NA,
                        group_a_mean=mean_a, group_b_mean=mean_b,
                        direction=NA, stringsAsFactors=FALSE)
    
    if (nrow(long) < 6) { rows_list[[i]] <- empty; next }
    
    tryCatch({
      fit    <- aov(value ~ condition * sex, data=long)
      aov_s  <- summary(fit)[[1]]
      rnames <- trimws(rownames(aov_s))
      cidx   <- which(rnames == "condition")
      omnibus_F <- if (length(cidx)) aov_s[["F value"]][cidx] else NA
      omnibus_p <- if (length(cidx)) aov_s[["Pr(>F)"]][cidx]  else NA
      
      em    <- emmeans::emmeans(fit, ~ condition)
      con   <- emmeans::contrast(em, method="pairwise", adjust="tukey")
      cdf   <- as.data.frame(con)
      ridx  <- pairwise_row(ga, gb, cond_names)   # select by position, not label
      crow  <- if (!is.na(ridx) && ridx <= nrow(cdf)) cdf[ridx, , drop=FALSE] else cdf[0, ]
      
      # Direction-safe estimate straight from the marginal means: B - A
      # (treatment - reference), independent of the contrast's internal ordering.
      emm_df <- as.data.frame(em)
      m_a    <- emm_df$emmean[emm_df$condition == ga]
      m_b    <- emm_df$emmean[emm_df$condition == gb]
      estimate <- if (length(m_a) && length(m_b)) m_b - m_a else NA
      se       <- if (nrow(crow)) crow$SE[1]      else NA
      p_tukey  <- if (nrow(crow)) crow$p.value[1] else NA
      t_s      <- if (!is.na(estimate) && !is.na(se) && se > 0) estimate / se else NA
      
      va <- c(as.numeric(norm_m[i, groups_m[[ga]], drop=TRUE]),
              as.numeric(norm_f[i, groups_f[[ga]], drop=TRUE]))
      vb <- c(as.numeric(norm_m[i, groups_m[[gb]], drop=TRUE]),
              as.numeric(norm_f[i, groups_f[[gb]], drop=TRUE]))
      va <- va[!is.na(va)]; vb <- vb[!is.na(vb)]
      n_a <- length(va); n_b <- length(vb)
      sd_p <- sqrt(((n_a-1)*var(va)+(n_b-1)*var(vb))/(n_a+n_b-2))
      cd   <- if (!is.na(sd_p) && sd_p > 0) (mean(vb)-mean(va))/sd_p else NA  # B - A
      fc   <- if (!is.na(mean_a) && mean_a > 0) mean_b/mean_a else NA
      dir  <- if (!is.na(mean_a) && !is.na(mean_b)) {
        if (mean_a > mean_b) paste(ga,">",gb) else paste(ga,"<",gb)
      } else NA
      
      rows_list[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                   comparison=paste(ga,"-",gb),
                                   omnibus_F=omnibus_F, omnibus_p=omnibus_p,
                                   contrast_F=if(!is.na(t_s)) t_s^2 else NA,
                                   estimate=estimate, std_error=se, t_statistic=t_s, p_tukey=p_tukey,
                                   fold_change=fc, cohens_d=cd,
                                   group_a_mean=mean_a, group_b_mean=mean_b,
                                   direction=dir, stringsAsFactors=FALSE)
    }, error=function(e) { rows_list[[i]] <<- empty })
  }
  
  res <- do.call(rbind, rows_list)
  valid_t <- !is.na(res$p_tukey);   res$p_fdr <- NA
  res$p_fdr[valid_t] <- p.adjust(res$p_tukey[valid_t], method="BH")
  valid_a <- !is.na(res$omnibus_p); res$omnibus_fdr <- NA
  res$omnibus_fdr[valid_a] <- p.adjust(res$omnibus_p[valid_a], method="BH")
  res$sig_tukey       <- sig_stars(res$p_tukey)
  res$sig_fdr         <- sig_stars(res$p_fdr)
  res$sig_omnibus_fdr <- sig_stars(res$omnibus_fdr)
  
  twoway_results[[label]] <- res
  cat(sprintf("    %-35s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}

twoway_results[["Sex Main Effect"]]          <- sex_effect_df
twoway_results[["Condition x Sex Interact"]] <- int_effect_df
write_excel(twoway_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_TwoWay.xlsx"))

# ==============================================================================
# SUPPLEMENTARY 2: Negative binomial - condition x sex, raw counts
# ==============================================================================
cat("\n--- SUPPLEMENTARY: Negative binomial (condition x sex, raw counts) ---\n")

negbin_results <- list()
for (pair in pairs) {
  ga <- pair[1]; gb <- pair[2]; label <- paste(ga, "vs", gb)
  cat(sprintf("  %s\n", label))
  
  rows_list <- vector("list", nrow(df_m))
  for (i in seq_len(nrow(df_m))) {
    acr <- df_m$acronym[i]; rgn <- df_m$name[i]
    vol <- df_m$volume_mm3[i]
    raw_a <- c(as.numeric(df_m[i, groups_m[[ga]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[ga]], drop=TRUE]))
    raw_b <- c(as.numeric(df_m[i, groups_m[[gb]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[gb]], drop=TRUE]))
    raw_a <- raw_a[!is.na(raw_a)]; raw_b <- raw_b[!is.na(raw_b)]
    mean_a <- if (length(raw_a)) mean(raw_a) else NA
    mean_b <- if (length(raw_b)) mean(raw_b) else NA
    
    empty <- data.frame(region_name=rgn, acronym=acr, comparison=paste(ga,"-",gb),
                        negbin_deviance=NA, negbin_p=NA, estimate=NA, std_error=NA,
                        z_statistic=NA, p_contrast=NA, fold_change=NA, IRR=NA,
                        group_a_mean=mean_a, group_b_mean=mean_b,
                        direction=NA, stringsAsFactors=FALSE)
    
    if (is.na(vol) || vol <= 0) { rows_list[[i]] <- empty; next }
    
    rows <- list()
    for (g in cond_names) {
      for (col in groups_m[[g]]) {
        cnt <- df_m[i, col]
        if (!is.na(cnt) && cnt >= 0)
          rows[[length(rows)+1]] <- data.frame(count=as.integer(round(cnt)),
                                               condition=g, sex="Male",   log_vol=log(vol), stringsAsFactors=FALSE)
      }
      for (col in groups_f[[g]]) {
        cnt <- df_f[i, col]
        if (!is.na(cnt) && cnt >= 0)
          rows[[length(rows)+1]] <- data.frame(count=as.integer(round(cnt)),
                                               condition=g, sex="Female", log_vol=log(vol), stringsAsFactors=FALSE)
      }
    }
    if (length(rows) < 8) { rows_list[[i]] <- empty; next }
    long <- do.call(rbind, rows)
    long$condition <- factor(long$condition, levels=cond_names)
    long$sex       <- factor(long$sex, levels=c("Male","Female"))
    
    tryCatch({
      fit      <- MASS::glm.nb(count ~ condition * sex + offset(log_vol), data=long)
      fit_null <- MASS::glm.nb(count ~ 1               + offset(log_vol), data=long)
      lrt      <- anova(fit_null, fit)
      negbin_p  <- lrt[["Pr(Chi)"]][2]
      negbin_dev <- lrt[["LR stat."]][2]
      
      em   <- emmeans::emmeans(fit, ~ condition)
      con  <- emmeans::contrast(em, method="pairwise", adjust="tukey")
      cdf  <- as.data.frame(con)
      ridx <- pairwise_row(ga, gb, cond_names)    # select by position, not label
      crow <- if (!is.na(ridx) && ridx <= nrow(cdf)) cdf[ridx, , drop=FALSE] else cdf[0, ]
      if (nrow(crow) == 0) stop("no matching contrast found")
      
      # Direction-safe log-rate difference from marginal means: B - A
      # (treatment - reference), so IRR = exp(est_log) is a treatment/reference
      # rate ratio and reads the same direction as fold_change.
      emm_df  <- as.data.frame(em)
      l_a     <- emm_df$emmean[emm_df$condition == ga]
      l_b     <- emm_df$emmean[emm_df$condition == gb]
      est_log <- if (length(l_a) && length(l_b)) l_b - l_a else NA
      se_log  <- crow$SE[1]
      z_stat  <- if (!is.na(est_log) && !is.na(se_log) && se_log > 0) est_log/se_log else NA
      p_con   <- crow$p.value[1]
      IRR     <- exp(est_log)
      fc  <- if (!is.na(mean_a) && mean_a > 0) mean_b/mean_a else NA
      dir <- if (!is.na(mean_a) && !is.na(mean_b)) {
        if (mean_a > mean_b) paste(ga,">",gb) else paste(ga,"<",gb)
      } else NA
      
      rows_list[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                   comparison=paste(ga,"-",gb), negbin_deviance=negbin_dev, negbin_p=negbin_p,
                                   estimate=est_log, std_error=se_log, z_statistic=z_stat, p_contrast=p_con,
                                   fold_change=fc, IRR=IRR, group_a_mean=mean_a, group_b_mean=mean_b,
                                   direction=dir, stringsAsFactors=FALSE)
    }, error=function(e) { rows_list[[i]] <<- empty })
  }
  
  res <- do.call(rbind, rows_list)
  valid_c <- !is.na(res$p_contrast); res$p_fdr <- NA
  res$p_fdr[valid_c] <- p.adjust(res$p_contrast[valid_c], method="BH")
  valid_m <- !is.na(res$negbin_p);   res$negbin_fdr <- NA
  res$negbin_fdr[valid_m] <- p.adjust(res$negbin_p[valid_m], method="BH")
  res$sig_tukey      <- sig_stars(res$p_contrast)
  res$sig_fdr        <- sig_stars(res$p_fdr)
  res$sig_negbin_fdr <- sig_stars(res$negbin_fdr)
  
  negbin_results[[label]] <- res
  cat(sprintf("    sig contrast: %3d (%3d)   FDR: %3d (%3d)\n",
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
write_excel(negbin_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_NegBin.xlsx"))

# ---- Save normalized data ----------------------------------------------------
write.csv(norm_comb, file.path(OUTPUT_DIR, "normalized_data_combined.csv"), row.names=FALSE)
write.csv(norm_m,    file.path(OUTPUT_DIR, "normalized_data_males.csv"),    row.names=FALSE)
write.csv(norm_f,    file.path(OUTPUT_DIR, "normalized_data_females.csv"),  row.names=FALSE)

# ---- Summary -----------------------------------------------------------------
cat("\n", strrep("=",65), "\n", sep="")
cat("SUMMARY\n")
cat(strrep("=",65), "\n", sep="")
cat("Normalization: log10(cells/mm3) | ANOVA on pooled males + females\n")
cat("Counts are regions significant at 0.05, with the trend-inclusive count\n")
cat("(0.10) in parentheses. Reported results use the 0.05 counts.\n\n")

cat("PRIMARY - One-way ANOVA (pooled):\n")
for (label in names(primary_results)) {
  res <- primary_results[[label]]
  cat(sprintf("  %-35s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
cat("\nSECONDARY - Sex-stratified ANOVA:\n")
for (label in names(bysex_results)) {
  res <- bysex_results[[label]]
  cat(sprintf("  %-40s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
cat("\nSUPPLEMENTARY - Two-way ANOVA (condition x sex):\n")
cat(sprintf("  Sex main effect:                    %3d (%3d) regions\n",
            n_below(sex_effect_df$q_sex, 0.05),
            n_below(sex_effect_df$q_sex, 0.10)))
cat(sprintf("  Condition x Sex interaction:        %3d (%3d) regions\n",
            n_below(int_effect_df$q_interaction, 0.05),
            n_below(int_effect_df$q_interaction, 0.10)))
for (label in setdiff(names(twoway_results),
                      c("Sex Main Effect","Condition x Sex Interact"))) {
  res <- twoway_results[[label]]
  cat(sprintf("  %-35s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
cat("\nOutputs:\n")
cat("  Drug_Statistical_Results_Primary.xlsx  [PRIMARY]\n")
cat("  Drug_Statistical_Results_BySex.xlsx    [SECONDARY]\n")
cat("  Drug_Statistical_Results_TwoWay.xlsx   [SUPPLEMENTARY]\n")
cat("  Drug_Statistical_Results_NegBin.xlsx   [SUPPLEMENTARY]\n")
cat("  normalized_data_combined/males/females.csv\n")
cat("\nDone.\n")