#!/usr/bin/env Rscript
# make_all_from_release.R  (PUBLIC, ENGLISH-ONLY, ONE-STOP)
#
# Purpose
#   Reproduce the main numeric tables and core figures from the *released TSVs*.
#   Designed for GitHub release packages (data/ + analysis/ + outputs/).
#
# Outputs (written to --outdir)
#   Tables:
#     Table2.tsv
#     Table3.tsv
#   Figures (separate PDFs, Illustrator-friendly):
#     Fig2A_CC_vs_lesionVol.pdf
#     Fig2B_wholebrain_mean_logJ.pdf
#     Fig3A_dDice_margins.pdf
#     Fig3B_dASSD_margins.pdf
#     Fig3C_lesion_mean_logJ.pdf
#     Fig4_dCC_vs_ventricle_fraction.pdf
#   Supplementary tables:
#     S4_subjectmask_fairness.tsv        (optional input)
#     S5_lesionVol_regression.tsv        (dCC, dlogJ_mean ~ lesion volume + covariates)
#     S6_absMeanLogJ_vs_ASSD.tsv         (mechanistic regression)
#     S7_intralesional_dlogJ_tertiles.tsv
#     S9_ventricle_fraction_dCC.tsv
#
# Inputs (default: <root>/data)
#   release_per_subject_metrics.tsv
#   lesion_roundtrip_dice_assd_native.tsv
#   lesion_logjac_volume_T1_MNI_EBT.tsv
# Optional:
#   bg_CC_gradCC_subjectmask_warp_table.tsv
#
# Notes
#   - Uses base pdf() (no cairo/X11 dependency)
#   - Uses ASCII-only plot labels to avoid Unicode rendering issues

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# -------------------- args --------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

ROOT <- get_arg("--root", "..")
DATA_DIR <- get_arg("--data-dir", file.path(ROOT, "data"))
OUTDIR <- get_arg("--outdir", file.path(ROOT, "outputs"))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

RELEASE_TSV <- get_arg("--release-tsv", file.path(DATA_DIR, "release_per_subject_metrics.tsv"))
LESION_RT_TSV <- get_arg("--lesion-rt-tsv", file.path(DATA_DIR, "lesion_roundtrip_dice_assd_native.tsv"))
LESION_LOGJ_TSV <- get_arg("--lesion-logj-tsv", file.path(DATA_DIR, "lesion_logjac_volume_T1_MNI_EBT.tsv"))
SUBMASK_SUMMARY <- get_arg("--submask-summary", file.path(DATA_DIR, "bg_CC_gradCC_subjectmask_warp_table.tsv"))

# -------- figure sizing / style knobs --------
FIG2A_W <- as.numeric(get_arg("--fig2a-w", "6.8"))
FIG2A_H <- as.numeric(get_arg("--fig2a-h", "4.5"))
FIG2B_W <- as.numeric(get_arg("--fig2b-w", "3.4"))
FIG2B_H <- as.numeric(get_arg("--fig2b-h", "4.5"))

FIG3_W  <- as.numeric(get_arg("--fig3-w",  "3.4"))
FIG3_H  <- as.numeric(get_arg("--fig3-h",  "4.5"))

FIG4_W  <- as.numeric(get_arg("--fig4-w",  "6.8"))
FIG4_H  <- as.numeric(get_arg("--fig4-h",  "4.8"))

BASE_SIZE <- as.numeric(get_arg("--base-size", "12"))
LINE_WIDTH <- as.numeric(get_arg("--line-width", "1.1"))

POINT_ALPHA <- as.numeric(get_arg("--point-alpha", "0.12"))
POINT_SIZE  <- as.numeric(get_arg("--point-size",  "1.5"))
POINT_STROKE<- as.numeric(get_arg("--point-stroke","0.4"))

JIT_ALPHA <- as.numeric(get_arg("--jit-alpha", "0.25"))
JIT_SIZE  <- as.numeric(get_arg("--jit-size",  "0.85"))
JIT_WIDTH <- as.numeric(get_arg("--jit-width", "0.18"))

BOX_LW <- as.numeric(get_arg("--box-lw", "0.6"))
BOX_W  <- as.numeric(get_arg("--box-w",  "0.6"))
BOX_FILL_STR <- get_arg("--box-fill", "NA") 
BOX_FILL <- if (toupper(BOX_FILL_STR) %in% c("NA","NONE","TRANSPARENT")) NA else BOX_FILL_STR

MARGIN_COL <- get_arg("--margin-col", "#d62728")
MARGIN_LW  <- as.numeric(get_arg("--margin-lw", "1.4"))
MARGIN_LTY <- get_arg("--margin-lty", "dotted")

COL_MNI <- get_arg("--col-mni", "#1f77b4")
COL_EBT <- get_arg("--col-ebt", "#d62728")

GRID_MAJOR_COL <- get_arg("--grid-major", "grey85")
GRID_MINOR_COL <- get_arg("--grid-minor", "grey92")
GRID_MAJOR_W   <- as.numeric(get_arg("--grid-major-w", "0.3"))
GRID_MINOR_W   <- as.numeric(get_arg("--grid-minor-w", "0.2"))

cat("[INFO] ROOT      :", ROOT, "\n")
cat("[INFO] DATA_DIR  :", DATA_DIR, "\n")
cat("[INFO] OUTDIR    :", OUTDIR, "\n\n")
cat("[INFO] RELEASE_TSV     :", RELEASE_TSV, "\n")
cat("[INFO] LESION_RT_TSV   :", LESION_RT_TSV, "\n")
cat("[INFO] LESION_LOGJ_TSV :", LESION_LOGJ_TSV, "\n")
cat("[INFO] SUBMASK_SUMMARY :", SUBMASK_SUMMARY, "\n\n")

if (!file.exists(RELEASE_TSV)) stop("[FATAL] release TSV not found: ", RELEASE_TSV)
if (!file.exists(LESION_RT_TSV)) stop("[FATAL] lesion round-trip TSV not found: ", LESION_RT_TSV)
if (!file.exists(LESION_LOGJ_TSV)) stop("[FATAL] lesion logJ TSV not found: ", LESION_LOGJ_TSV)

# -------------------- helpers --------------------
require_cols <- function(df, cols, name) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop("[FATAL] missing columns in ", name, ": ", paste(miss, collapse = ", "))
}

fmt_mean_sd <- function(m, s, digits = 3) {
  sprintf(paste0("%.", digits, "f\u00B1%.", digits, "f"), m, s)
}
fmt_ci <- function(l, u, digits = 3) {
  paste0(sprintf(paste0("%.", digits, "f"), l), " to ", sprintf(paste0("%.", digits, "f"), u))
}
fmt_delta_ci <- function(delta, l, u, digits = 3) {
  paste0(sprintf(paste0("%.", digits, "f"), delta), " (", fmt_ci(l, u, digits), ")")
}
fmt_p <- function(p) {
  vapply(p, function(pp) {
    if (is.na(pp)) return(NA_character_)
    if (pp < 0.001) return("<0.001")
    sprintf("%.3f", pp)
  }, character(1))
}

# Robust numeric parser for voxel size strings (e.g., "1", "1.0", "1x1x1", "1.0 1.0 1.0")
parse_resolution <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  s <- as.character(x)
  s <- gsub(",", ".", s)
  m <- regmatches(s, gregexpr("[0-9]+\\.?[0-9]*", s))
  out <- sapply(m, function(v) if (length(v) == 0) NA_real_ else as.numeric(v[1]))
  as.numeric(out)
}

theme_masume <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(colour = GRID_MAJOR_COL, linewidth = GRID_MAJOR_W),
      panel.grid.minor = element_line(colour = GRID_MINOR_COL, linewidth = GRID_MINOR_W),
      panel.border = element_rect(colour = "black", linewidth = 0.6),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(6, 6, 6, 6)
    )
}

safe_pdf <- function(path, w, h) {
  pdf(path, width = w, height = h, useDingbats = FALSE)
}

# -------------------- load data --------------------
rel <- read_tsv(RELEASE_TSV, show_col_types = FALSE)
lrt <- read_tsv(LESION_RT_TSV, show_col_types = FALSE)
llj <- read_tsv(LESION_LOGJ_TSV, show_col_types = FALSE)

# Required columns in release
require_cols(rel, c(
  "subject",
  "CC_MNI_bg", "CC_EBT_bg", "gCC_MNI_bg", "gCC_EBT_bg",
  "logJac_mean_MNI", "logJac_mean_EBT",
  "logJac_median_MNI", "logJac_median_EBT",
  "logJac_p5_MNI", "logJac_p5_EBT",
  "logJac_p95_MNI", "logJac_p95_EBT",
  "Vlesion_native_mm3",
  "ventricle_fraction",
  "ScannerBrand", "ImageResolution"
), "release_per_subject_metrics")

# Required columns in lesion round-trip
# (ASSD columns are lower-case in the current public TSV)
require_cols(lrt, c("subject", "dice_MNI", "dice_EBT", "assd_MNI_mm", "assd_EBT_mm"), "lesion_roundtrip_dice_assd_native")

# Required columns in lesion logJ
require_cols(llj, c(
  "subject",
  "logJac_mean_MNI_lesion", "logJac_mean_EBT_lesion",
  "lesion_vol_T1_mm3"
), "lesion_logjac_volume_T1_MNI_EBT")

# -------------------- Table 2 --------------------
metric_specs <- tribble(
  ~Metric, ~mni_col, ~ebt_col, ~digits,
  "Background CC",        "CC_MNI_bg",        "CC_EBT_bg",        3,
  "Gradient-CC",          "gCC_MNI_bg",       "gCC_EBT_bg",       3,
  "log-Jacobian mean",    "logJac_mean_MNI",  "logJac_mean_EBT",  3,
  "log-Jacobian median",  "logJac_median_MNI","logJac_median_EBT",3,
  "log-Jacobian p5",      "logJac_p5_MNI",    "logJac_p5_EBT",    3,
  "log-Jacobian p95",     "logJac_p95_MNI",   "logJac_p95_EBT",   3
)

calc_pair_row <- function(df, mni_col, ebt_col, digits = 3) {
  x <- as.numeric(df[[mni_col]])
  y <- as.numeric(df[[ebt_col]])
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  d <- y - x
  tt <- t.test(y, x, paired = TRUE)
  list(
    n = length(d),
    mni_mean = mean(x), mni_sd = sd(x),
    ebt_mean = mean(y), ebt_sd = sd(y),
    d_mean = mean(d),
    ci_l = unname(tt$conf.int[1]),
    ci_u = unname(tt$conf.int[2]),
    p = unname(tt$p.value)
  )
}

rows2 <- metric_specs %>%
  rowwise() %>%
  mutate(tmp = list(calc_pair_row(rel, mni_col, ebt_col, digits))) %>%
  ungroup()

p_raw <- sapply(rows2$tmp, function(z) z$p)
p_holm <- p.adjust(p_raw, method = "holm")

Table2 <- rows2 %>%
  mutate(
    n = sapply(tmp, `[[`, "n"),
    `MNI152` = mapply(function(z, dig) fmt_mean_sd(z$mni_mean, z$mni_sd, dig), tmp, digits),
    `EBT`    = mapply(function(z, dig) fmt_mean_sd(z$ebt_mean, z$ebt_sd, dig), tmp, digits),
    `Delta (95% CI)` = mapply(function(z, dig) fmt_delta_ci(z$d_mean, z$ci_l, z$ci_u, dig), tmp, digits),
    `p-value` = sapply(p_holm, fmt_p)
  ) %>%
  select(Metric, MNI152, EBT, `Delta (95% CI)`, `p-value`)

write_tsv(Table2, file.path(OUTDIR, "Table2.tsv"))
cat("[INFO] wrote:", file.path(OUTDIR, "Table2.tsv"), "\n")

# -------------------- Table 3 --------------------
# Dice / ASSD absolute + delta + p + equivalence (TOST via 90% CI inclusion)
calc_equiv <- function(d, margin) {
  d <- d[is.finite(d)]
  n <- length(d)
  m <- mean(d)
  s <- sd(d)
  se <- s / sqrt(n)
  tcrit <- qt(0.95, df = n - 1)  # 90% CI
  ci90_l <- m - tcrit * se
  ci90_u <- m + tcrit * se
  equiv <- (ci90_l >= -margin) && (ci90_u <= margin)
  list(ci90_l = ci90_l, ci90_u = ci90_u, equiv = equiv)
}

dice_mni <- as.numeric(lrt$dice_MNI)
dice_ebt <- as.numeric(lrt$dice_EBT)
assd_mni <- as.numeric(lrt$assd_MNI_mm)
assd_ebt <- as.numeric(lrt$assd_EBT_mm)

ok_dice <- is.finite(dice_mni) & is.finite(dice_ebt)
ok_assd <- is.finite(assd_mni) & is.finite(assd_ebt)

d_dice <- dice_ebt[ok_dice] - dice_mni[ok_dice]
d_assd <- assd_ebt[ok_assd] - assd_mni[ok_assd]

t_dice <- t.test(dice_ebt[ok_dice], dice_mni[ok_dice], paired = TRUE)
t_assd <- t.test(assd_ebt[ok_assd], assd_mni[ok_assd], paired = TRUE)

# intralesional mean logJ difference (EBT - MNI)
llj2 <- llj %>%
  transmute(
    subject,
    meanLogJ_MNI_les = as.numeric(logJac_mean_MNI_lesion),
    meanLogJ_EBT_les = as.numeric(logJac_mean_EBT_lesion)
  ) %>%
  filter(is.finite(meanLogJ_MNI_les), is.finite(meanLogJ_EBT_les))

d_logj_les <- llj2$meanLogJ_EBT_les - llj2$meanLogJ_MNI_les

t_logj_les <- t.test(d_logj_les, mu = 0)

# equivalence checks
EQ_DICE <- calc_equiv(d_dice, margin = 0.02)
EQ_ASSD <- calc_equiv(d_assd, margin = 0.10)

Table3 <- tibble(
  Measure = c("Dice", "ASSD (mm)", "Lesion logJ mean"),
  `MNI152 (mean±SD)` = c(
    fmt_mean_sd(mean(dice_mni[ok_dice]), sd(dice_mni[ok_dice]), 3),
    fmt_mean_sd(mean(assd_mni[ok_assd]), sd(assd_mni[ok_assd]), 3),
    "—"
  ),
  `EBT (mean±SD)` = c(
    fmt_mean_sd(mean(dice_ebt[ok_dice]), sd(dice_ebt[ok_dice]), 3),
    fmt_mean_sd(mean(assd_ebt[ok_assd]), sd(assd_ebt[ok_assd]), 3),
    "—"
  ),
  `Delta (95% CI)` = c(
    fmt_delta_ci(mean(d_dice), unname(t_dice$conf.int[1]), unname(t_dice$conf.int[2]), 3),
    fmt_delta_ci(mean(d_assd), unname(t_assd$conf.int[1]), unname(t_assd$conf.int[2]), 3),
    fmt_delta_ci(mean(d_logj_les), unname(t_logj_les$conf.int[1]), unname(t_logj_les$conf.int[2]), 3)
  ),
  `p-value` = c(fmt_p(t_dice$p.value), fmt_p(t_assd$p.value), fmt_p(t_logj_les$p.value)),
  `Equiv. (TOST)` = c(
    ifelse(EQ_DICE$equiv, "YES (±0.02)", "NO (±0.02)"),
    ifelse(EQ_ASSD$equiv, "YES (±0.10)", "NO (±0.10)"),
    "N/A"
  )
)

write_tsv(Table3, file.path(OUTDIR, "Table3.tsv"))
cat("[INFO] wrote:", file.path(OUTDIR, "Table3.tsv"), "\n")

# -------------------- Figure 2A: CC vs lesion volume --------------------
# Build long CC
ccA <- rel %>%
  transmute(
    subject,
    log_lesion = log1p(as.numeric(Vlesion_native_mm3)),
    CC_MNI = as.numeric(CC_MNI_bg),
    CC_EBT = as.numeric(CC_EBT_bg)
  ) %>%
  filter(is.finite(log_lesion), is.finite(CC_MNI), is.finite(CC_EBT))

ccA_long <- ccA %>%
  pivot_longer(cols = c(CC_MNI, CC_EBT), names_to = "Template", values_to = "CC") %>%
  mutate(Template = ifelse(Template == "CC_MNI", "MNI152", "EBT"))

p2a <- ggplot(ccA_long, aes(x = log_lesion, y = CC, colour = Template, fill = Template)) +
  geom_point(shape = 21, alpha = POINT_ALPHA, size = POINT_SIZE, stroke = POINT_STROKE) +
  geom_smooth(method = "lm", se = FALSE, linewidth = LINE_WIDTH) +
  labs(x = "log(1 + lesion volume [mm^3])", y = "CC") +
  scale_colour_manual(values = c("MNI152" = COL_MNI, "EBT" = COL_EBT)) +
  scale_fill_manual(values   = c("MNI152" = COL_MNI, "EBT" = COL_EBT)) +
  theme_masume(BASE_SIZE) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0.0, 1.0, by = 0.1), labels = function(x) sprintf("%.1f", x)) +
  scale_x_continuous(breaks = seq(4, 12, by = 2))

out2a <- file.path(OUTDIR, "Fig2A_CC_vs_lesionVol.pdf")
safe_pdf(out2a, FIG2A_W, FIG2A_H)
print(p2a)
dev.off()
cat("[INFO] wrote:", out2a, "\n")

# -------------------- Figure 2B: whole-brain mean logJ --------------------
ljB <- rel %>%
  transmute(
    subject,
    logJ_MNI = as.numeric(logJac_mean_MNI),
    logJ_EBT = as.numeric(logJac_mean_EBT)
  ) %>%
  filter(is.finite(logJ_MNI), is.finite(logJ_EBT))

ljB_long <- ljB %>%
  pivot_longer(cols = c(logJ_MNI, logJ_EBT), names_to = "Template", values_to = "logJ") %>%
  mutate(Template = ifelse(Template == "logJ_MNI", "MNI152", "EBT"))

p2b <- ggplot(ljB_long, aes(x = Template, y = logJ, colour = Template)) +
  geom_jitter(alpha = JIT_ALPHA, width = JIT_WIDTH, size = JIT_SIZE) +
  geom_boxplot(outlier.shape = NA, linewidth = BOX_LW, colour = "black", width = BOX_W, fill = BOX_FILL) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  labs(x = NULL, y = "Whole-brain mean logJ") +
  scale_colour_manual(values = c("MNI152" = COL_MNI, "EBT" = COL_EBT)) +
  theme_masume(BASE_SIZE) +
  theme(legend.position = "none")

out2b <- file.path(OUTDIR, "Fig2B_wholebrain_mean_logJ.pdf")
safe_pdf(out2b, FIG2B_W, FIG2B_H)
print(p2b)
dev.off()
cat("[INFO] wrote:", out2b, "\n")

# -------------------- Figure 3A/B/C: lesion round-trip + intralesional logJ --------------------
# Robust column picker (case-insensitive)
pick_col <- function(df, candidates) {
  nms <- names(df)
  nms_l <- tolower(nms)
  for (cand in candidates) {
    idx <- which(nms_l == tolower(cand))
    if (length(idx) == 1) return(nms[idx])
  }
  return(NULL)
}
require_cols2 <- function(df, cols, label) {
  miss <- cols[is.na(cols) | !(cols %in% names(df))]
  if (length(miss) > 0) stop(sprintf("[FATAL] missing columns in %s: %s", label, paste(miss, collapse=", ")))
}

# ----- style defaults if not defined earlier -----
if (!exists("JIT_ALPHA"))  JIT_ALPHA <- 0.25
if (!exists("JIT_SIZE"))   JIT_SIZE  <- 0.85
if (!exists("JIT_WIDTH"))  JIT_WIDTH <- 0.25
if (!exists("BOX_LW"))     BOX_LW    <- 0.6
if (!exists("BOX_W"))      BOX_W     <- 0.35
if (!exists("BOX_FILL"))   BOX_FILL  <- NA
if (!exists("BASE_SIZE"))  BASE_SIZE <- 12
MARGIN_COL <- "#d62728"
MARGIN_LW  <- 0.6
MARGIN_LTY <- "dashed"

# Figure size
FIG3_W <- as.numeric(get_arg("--fig3-w", "3.4"))
FIG3_H <- as.numeric(get_arg("--fig3-h", "4.5"))

# Axis controls (old-style look)
FIG3A_YMIN <- as.numeric(get_arg("--fig3a-ymin", "-0.065"))
FIG3A_YMAX <- as.numeric(get_arg("--fig3a-ymax",  "0.065"))
FIG3B_YMIN <- as.numeric(get_arg("--fig3b-ymin", "-0.20"))
FIG3B_YMAX <- as.numeric(get_arg("--fig3b-ymax",  "0.20"))
FIG3_JIT_WIDTH <- as.numeric(get_arg("--fig3-jit-width", "0.35"))

# -------------------- Fig3A / Fig3B from lesion round-trip TSV --------------------
# Detect columns in lesion_roundtrip_dice_assd_native.tsv (lrt)
lrt_id_col   <- pick_col(lrt, c("subject", "case", "id"))
dice_mni_col <- pick_col(lrt, c("dice_mni", "dice_mni152", "dice_MNI", "dice_MNI152"))
dice_ebt_col <- pick_col(lrt, c("dice_ebt", "dice_EBT"))
assd_mni_col <- pick_col(lrt, c("assd_mni_mm", "assd_mni", "ASSD_MNI_mm", "ASSD_MNI"))
assd_ebt_col <- pick_col(lrt, c("assd_ebt_mm", "assd_ebt", "ASSD_EBT_mm", "ASSD_EBT"))

require_cols2(lrt,
              c(lrt_id_col, dice_mni_col, dice_ebt_col, assd_mni_col, assd_ebt_col),
              "lesion_roundtrip_dice_assd_native")

fig3_df <- lrt %>%
  transmute(
    subject = .data[[lrt_id_col]],
    dDice   = as.numeric(.data[[dice_ebt_col]]) - as.numeric(.data[[dice_mni_col]]),
    dASSD   = as.numeric(.data[[assd_ebt_col]]) - as.numeric(.data[[assd_mni_col]])
  ) %>%
  filter(is.finite(dDice), is.finite(dASSD))

if (nrow(fig3_df) < 10) stop("[FATAL] Too few rows for Fig3A/B")

# ---- Fig3A: delta Dice ----
p3a <- ggplot(fig3_df, aes(x = 1, y = dDice)) +
  geom_jitter(width = FIG3_JIT_WIDTH, height = 0,
              alpha = JIT_ALPHA, size = JIT_SIZE, colour = "grey40") +
  geom_boxplot(width = BOX_W, outlier.shape = NA,
               linewidth = BOX_LW, fill = BOX_FILL, colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5, colour = "black") +
  geom_hline(yintercept = c(-0.02, 0.02),
             linetype = MARGIN_LTY, linewidth = MARGIN_LW, colour = MARGIN_COL) +
  labs(x = NULL, y = "Delta Dice (EBT - MNI152)") +
  theme_bw(base_size = BASE_SIZE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(
    breaks = c(-0.06, -0.03, 0.00, 0.03, 0.06),
    labels = function(x) sprintf("%.2f", x)
  ) +
  coord_cartesian(ylim = c(FIG3A_YMIN, FIG3A_YMAX))

out3a <- file.path(OUTDIR, "Fig3A_dDice_margins.pdf")
safe_pdf(out3a, FIG3_W, FIG3_H); print(p3a); dev.off()
cat("[INFO] wrote:", out3a, "\n")

# ---- Fig3B: delta ASSD ----
p3b <- ggplot(fig3_df, aes(x = 1, y = dASSD)) +
  geom_jitter(width = FIG3_JIT_WIDTH, height = 0,
              alpha = JIT_ALPHA, size = JIT_SIZE, colour = "grey40") +
  geom_boxplot(width = BOX_W, outlier.shape = NA,
               linewidth = BOX_LW, fill = BOX_FILL, colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5, colour = "black") +
  geom_hline(yintercept = c(-0.10, 0.10),
             linetype = MARGIN_LTY, linewidth = MARGIN_LW, colour = MARGIN_COL) +
  labs(x = NULL, y = "Delta ASSD (mm; EBT - MNI152)") +
  theme_bw(base_size = BASE_SIZE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(
    breaks = c(-0.2,-0.1, 0.0, 0.1, 0.2),
    labels = function(x) sprintf("%.1f", x)
  ) +
  coord_cartesian(ylim = c(FIG3B_YMIN, FIG3B_YMAX))

out3b <- file.path(OUTDIR, "Fig3B_dASSD_margins.pdf")
safe_pdf(out3b, FIG3_W, FIG3_H); print(p3b); dev.off()
cat("[INFO] wrote:", out3b, "\n")

# -------------------- Fig3C: lesion mean logJ (absolute per-template values) --------------------
# Detect columns in lesion_logjac_volume_T1_MNI_EBT.tsv (llj)
llj_id_col <- pick_col(llj, c("subject", "case", "id"))
mni_lcol   <- pick_col(llj, c("logjac_mean_mni_lesion", "logjac_mean_mni", "logJac_mean_MNI_lesion"))
ebt_lcol   <- pick_col(llj, c("logjac_mean_ebt_lesion", "logjac_mean_ebt", "logJac_mean_EBT_lesion"))

require_cols2(llj, c(llj_id_col, mni_lcol, ebt_lcol), "lesion_logjac_volume_T1_MNI_EBT")

fig3c_df <- llj %>%
  transmute(
    subject = .data[[llj_id_col]],
    MNI152  = as.numeric(.data[[mni_lcol]]),
    EBT     = as.numeric(.data[[ebt_lcol]])
  ) %>%
  filter(is.finite(MNI152), is.finite(EBT))

if (nrow(fig3c_df) < 10) stop("[FATAL] Too few rows for Fig3C")

fig3c_long <- fig3c_df %>%
  pivot_longer(cols = c(MNI152, EBT), names_to = "Template", values_to = "logJ") %>%
  mutate(Template = factor(Template, levels = c("MNI152", "EBT")))

FIG3C_W <- as.numeric(get_arg("--fig3c-w", "4.5"))
FIG3C_H <- as.numeric(get_arg("--fig3c-h", "4.5"))

# Colors (use same global if you have; otherwise default)
if (!exists("COL_MNI")) COL_MNI <- "#1f77b4"
if (!exists("COL_EBT")) COL_EBT <- "#d62728"

p3c <- ggplot(fig3c_long, aes(x = Template, y = logJ, colour = Template)) +
  geom_jitter(alpha = JIT_ALPHA, width = JIT_WIDTH, size = JIT_SIZE) +
  geom_boxplot(outlier.shape = NA, linewidth = BOX_LW, colour = "black",
               width = 0.6, fill = BOX_FILL) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = c("MNI152" = COL_MNI, "EBT" = COL_EBT)) +
  labs(x = NULL, y = "Lesion mean logJ") +
  theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = "none")

out3c <- file.path(OUTDIR, "Fig3C_lesion_mean_logJ.pdf")
safe_pdf(out3c, FIG3C_W, FIG3C_H); print(p3c); dev.off()
cat("[INFO] wrote:", out3c, "\n")


# -------------------- Figure 4: dCC vs ventricle fraction --------------------
fig4_df <- rel %>%
  transmute(
    subject,
    ventfrac = as.numeric(ventricle_fraction),
    dCC = as.numeric(CC_EBT_bg) - as.numeric(CC_MNI_bg)
  ) %>%
  filter(is.finite(ventfrac), is.finite(dCC))

if (nrow(fig4_df) < 10) stop("[FATAL] Too few rows for Figure 4")

FIG4_W <- as.numeric(get_arg("--fig4-w", "6.8"))
FIG4_H <- as.numeric(get_arg("--fig4-h", "4.5"))

p4 <- ggplot(fig4_df, aes(x = ventfrac, y = dCC)) +
  geom_point(alpha = 0.25, size = 1.2, colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, linewidth = LINE_WIDTH, colour = "black") +
  labs(
    x = "Ventricle fraction (ventricular volume / ICV)",
    y = "dCC (EBT - MNI152)"
  ) +
  theme_bw(base_size = BASE_SIZE)

out4 <- file.path(OUTDIR, "Fig4_dCC_vs_ventricle_fraction.pdf")
safe_pdf(out4, FIG4_W, FIG4_H)
print(p4)
dev.off()
cat("[INFO] wrote:", out4, "\n")

# -------------------- Supplementary Table S4 (subject-mask fairness) --------------------
if (file.exists(SUBMASK_SUMMARY)) {
  s4 <- read_tsv(SUBMASK_SUMMARY, show_col_types = FALSE)
  write_tsv(s4, file.path(OUTDIR, "S4_subjectmask_fairness.tsv"))
  cat("[INFO] wrote:", file.path(OUTDIR, "S4_subjectmask_fairness.tsv"), "\n")
} else {
  cat("[INFO] S4 skipped (file not found):", SUBMASK_SUMMARY, "\n")
}

# -------------------- Supplementary Table S5 (lesion volume + covariates) --------------------
# outcomes: dCC and dlogJ_mean; coefficient reported is the lesion-volume term
s5_df <- rel %>%
  transmute(
    subject,
    lesion_log = log1p(as.numeric(Vlesion_native_mm3)),
    ScannerBrand = as.factor(ScannerBrand),
    ImageResolution = parse_resolution(ImageResolution),
    dCC = as.numeric(CC_EBT_bg) - as.numeric(CC_MNI_bg),
    dlogJ_mean = as.numeric(logJac_mean_EBT) - as.numeric(logJac_mean_MNI)
  ) %>%
  filter(is.finite(lesion_log), is.finite(ImageResolution), is.finite(dCC), is.finite(dlogJ_mean))

build_formula <- function(yname, df) {
  terms <- c("lesion_log")
  if (nlevels(df$ScannerBrand) >= 2) terms <- c(terms, "ScannerBrand")
  if (length(unique(df$ImageResolution[is.finite(df$ImageResolution)])) >= 2) terms <- c(terms, "ImageResolution")
  as.formula(paste0(yname, " ~ ", paste(terms, collapse = " + ")))
}

fit_s5 <- function(yname) {
  f <- build_formula(yname, s5_df)
  m <- lm(f, data = s5_df)
  co <- summary(m)$coefficients
  # lesion_log row always exists
  beta <- unname(co["lesion_log", "Estimate"])
  se   <- unname(co["lesion_log", "Std. Error"])
  pval <- unname(co["lesion_log", "Pr(>|t|)"])
  dfree <- df.residual(m)
  tcrit <- qt(0.975, df = dfree)
  ci_l <- beta - tcrit * se
  ci_u <- beta + tcrit * se
  adjr2 <- unname(summary(m)$adj.r.squared)
  tibble(
    Outcome = yname,
    `Beta (per 1 log unit, 95% CI)` = fmt_delta_ci(beta, ci_l, ci_u, 3),
    `p-value` = fmt_p(pval),
    `Adjusted R2` = sprintf("%.3f", adjr2)
  )
}

S5 <- bind_rows(fit_s5("dCC"), fit_s5("dlogJ_mean"))
write_tsv(S5, file.path(OUTDIR, "S5_lesionVol_regression.tsv"))
cat("[INFO] wrote:", file.path(OUTDIR, "S5_lesionVol_regression.tsv"), "\n")

# -------------------- Supplementary Table S6 (abs mean logJ vs ASSD) --------------------
# MNI / EBT:
#   predictor = |mean logJ within lesion|  (template-specific)
#   outcome   = ASSD (template-specific round-trip)
# MNI152 - EBT:
#   predictor = dAbsMeanLogJ = |mean logJ|_MNI - |mean logJ|_EBT
#   outcome   = dASSD        = ASSD_MNI - ASSD_EBT
#
# Notes:
# - absMeanLogJ_* is |mean logJ| (absolute value of the intralesional mean logJ).
# - The difference row uses differences of the absolute summaries (NOT abs of the difference).
#   This keeps the "absolute" definition consistent across rows.

s6_df <- lrt %>%
  transmute(
    subject,
    ASSD_MNI = as.numeric(assd_MNI_mm),
    ASSD_EBT = as.numeric(assd_EBT_mm)
  ) %>%
  inner_join(llj2, by = "subject") %>%
  mutate(
    absMeanLogJ_MNI = abs(as.numeric(meanLogJ_MNI_les)),
    absMeanLogJ_EBT = abs(as.numeric(meanLogJ_EBT_les)),
    dAbsMeanLogJ    = absMeanLogJ_MNI - absMeanLogJ_EBT,
    dASSD           = ASSD_MNI - ASSD_EBT
  ) %>%
  filter(
    is.finite(ASSD_MNI), is.finite(ASSD_EBT),
    is.finite(absMeanLogJ_MNI), is.finite(absMeanLogJ_EBT),
    is.finite(dAbsMeanLogJ), is.finite(dASSD)
  )

fit_line <- function(x, y) {
  m <- lm(y ~ x)
  co <- summary(m)$coefficients
  tibble(
    intercept = unname(co["(Intercept)", "Estimate"]),
    slope     = unname(co["x", "Estimate"]),
    R2        = unname(summary(m)$r.squared),
    p_slope   = unname(co["x", "Pr(>|t|)"])
  )
}

S6 <- bind_rows(
  bind_cols(
    tibble(model = "MNI", n = nrow(s6_df)),
    fit_line(s6_df$absMeanLogJ_MNI, s6_df$ASSD_MNI),
    tibble(
      spearman_rho = unname(cor(s6_df$absMeanLogJ_MNI, s6_df$ASSD_MNI, method = "spearman")),
      p_spearman   = unname(cor.test(s6_df$absMeanLogJ_MNI, s6_df$ASSD_MNI, method = "spearman", exact = FALSE)$p.value)
    )
  ),
  bind_cols(
    tibble(model = "EBT", n = nrow(s6_df)),
    fit_line(s6_df$absMeanLogJ_EBT, s6_df$ASSD_EBT),
    tibble(
      spearman_rho = unname(cor(s6_df$absMeanLogJ_EBT, s6_df$ASSD_EBT, method = "spearman")),
      p_spearman   = unname(cor.test(s6_df$absMeanLogJ_EBT, s6_df$ASSD_EBT, method = "spearman", exact = FALSE)$p.value)
    )
  ),
  bind_cols(
    tibble(model = "MNI152 - EBT", n = nrow(s6_df)),
    fit_line(s6_df$dAbsMeanLogJ, s6_df$dASSD),
    tibble(
      spearman_rho = unname(cor(s6_df$dAbsMeanLogJ, s6_df$dASSD, method = "spearman")),
      p_spearman   = unname(cor.test(s6_df$dAbsMeanLogJ, s6_df$dASSD, method = "spearman", exact = FALSE)$p.value)
    )
  )
) %>%
  mutate(
    intercept    = sprintf("%.4f", as.numeric(intercept)),
    slope        = sprintf("%.4f", as.numeric(slope)),
    R2           = sprintf("%.3f",  as.numeric(R2)),
    p_slope      = fmt_p(as.numeric(p_slope)),
    spearman_rho = sprintf("%.3f",  as.numeric(spearman_rho)),
    p_spearman   = fmt_p(as.numeric(p_spearman))
  )

write_tsv(S6, file.path(OUTDIR, "S6_absMeanLogJ_vs_ASSD.tsv"))
cat("[INFO] wrote:", file.path(OUTDIR, "S6_absMeanLogJ_vs_ASSD.tsv"), "\n")

# quick sanity checks (optional; safe to keep)
cat("[CHECK] absMeanLogJ_MNI min/max:", min(s6_df$absMeanLogJ_MNI), max(s6_df$absMeanLogJ_MNI), "\n")
cat("[CHECK] absMeanLogJ_EBT min/max:", min(s6_df$absMeanLogJ_EBT), max(s6_df$absMeanLogJ_EBT), "\n")
cat("[CHECK] any negative absMeanLogJ_MNI?", any(s6_df$absMeanLogJ_MNI < 0), "\n")
cat("[CHECK] any negative absMeanLogJ_EBT?", any(s6_df$absMeanLogJ_EBT < 0), "\n")

# -------------------- Supplementary Table S7 (intralesional dlogJ by lesion-volume tertiles) --------------------
# Use lesion_vol_T1_mm3 tertiles; outcome is (EBT - MNI) intralesional mean logJ
s7_df <- llj %>%
  transmute(
    subject,
    lesion_vol = as.numeric(lesion_vol_T1_mm3),
    dlogJ_les = as.numeric(logJac_mean_EBT_lesion) - as.numeric(logJac_mean_MNI_lesion)
  ) %>%
  filter(is.finite(lesion_vol), is.finite(dlogJ_les))

q <- quantile(s7_df$lesion_vol, probs = c(1/3, 2/3), na.rm = TRUE)

s7_df <- s7_df %>%
  mutate(group = case_when(
    lesion_vol <= q[1] ~ "small",
    lesion_vol <= q[2] ~ "medium",
    TRUE ~ "large"
  ))

one_group <- function(df) {
  n <- nrow(df)
  m <- mean(df$dlogJ_les)
  s <- sd(df$dlogJ_les)
  tt <- t.test(df$dlogJ_les, mu = 0)
  tibble(
    n = n,
    `dlogJ mean (mean±SD)` = fmt_mean_sd(m, s, 4),
    `dlogJ mean (95% CI)` = fmt_ci(unname(tt$conf.int[1]), unname(tt$conf.int[2]), 4),
    p_raw = unname(tt$p.value)
  )
}

S7_tmp <- s7_df %>%
  group_by(group) %>%
  group_modify(~one_group(.x)) %>%
  ungroup()

p_h <- p.adjust(S7_tmp$p_raw, method = "holm")
S7 <- S7_tmp %>%
  mutate(`p (Holm)` = sapply(p_h, fmt_p)) %>%
  select(group, n, `dlogJ mean (mean±SD)`, `dlogJ mean (95% CI)`, `p (Holm)`)

write_tsv(S7, file.path(OUTDIR, "S7_intralesional_dlogJ_tertiles.tsv"))
cat("[INFO] wrote:", file.path(OUTDIR, "S7_intralesional_dlogJ_tertiles.tsv"), "\n")

# -------------------- Supplementary Table S9 (ventricle fraction vs dCC) --------------------
# Simple regression: dCC ~ ventricle_fraction
s9_df <- rel %>%
  transmute(
    ventfrac = as.numeric(ventricle_fraction),
    dCC = as.numeric(CC_EBT_bg) - as.numeric(CC_MNI_bg)
  ) %>%
  filter(is.finite(ventfrac), is.finite(dCC))

m9 <- lm(dCC ~ ventfrac, data = s9_df)
co9 <- summary(m9)$coefficients
beta <- unname(co9["ventfrac", "Estimate"])
se   <- unname(co9["ventfrac", "Std. Error"])
pval <- unname(co9["ventfrac", "Pr(>|t|)"])

tcrit <- qt(0.975, df = df.residual(m9))
ci_l <- beta - tcrit * se
ci_u <- beta + tcrit * se

rho <- unname(cor(s9_df$ventfrac, s9_df$dCC, method = "spearman"))
rho_p <- unname(cor.test(s9_df$ventfrac, s9_df$dCC, method = "spearman", exact = FALSE)$p.value)

S9 <- tibble(
  n = nrow(s9_df),
  Outcome = "dCC (CC_EBT - CC_MNI152)",
  Predictor = "Ventricle fraction",
  `beta (95% CI)` = fmt_delta_ci(beta, ci_l, ci_u, 2),
  `Adjusted R2` = sprintf("%.3f", unname(summary(m9)$adj.r.squared)),
  `p-value` = fmt_p(pval),
  `Spearman rho` = sprintf("%.3f", rho),
  `p (Spearman)` = fmt_p(rho_p)
)

write_tsv(S9, file.path(OUTDIR, "S9_ventricle_fraction_dCC.tsv"))
cat("[INFO] wrote:", file.path(OUTDIR, "S9_ventricle_fraction_dCC.tsv"), "\n")

cat("[INFO] Done. Wrote outputs to:", OUTDIR, "\n")
