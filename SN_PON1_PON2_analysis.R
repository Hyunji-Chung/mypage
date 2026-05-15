# =============================================================================
# Dataset : SN protein data.xlsx  (pon2 branch)
# Study   : Parkinson's Disease — Substantia Nigra Proteomics
# Method  : Normalized intensity | Healthy vs PD
#
# Analysis:
#   1. limma DE (PD vs Healthy)
#   2. PON1 + PON2 boxplot (side-by-side) with p-value annotation
#   3. Full volcano plot with PON1 & PON2 highlighted
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", ask = FALSE)
  for (pkg in c("readxl", "tidyverse", "ggrepel", "ggpubr", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg)
  }
  library(limma)
  library(readxl)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(patchwork)
})

# ── 1. File paths ─────────────────────────────────────────────────────────────
INPUT_FILE <- "SN protein data.xlsx"
OUTPUT_DIR <- "SN_PON_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== PD SN Proteomics | PON1 & PON2 Expression Analysis ===\n")
cat("Input:", INPUT_FILE, "\n\n")

# ── 2. Sheet detection ────────────────────────────────────────────────────────
all_sheets <- excel_sheets(INPUT_FILE)
cat("Available sheets:", paste(all_sheets, collapse = ", "), "\n")

# Prefer "Normalized" sheet; fall back to first sheet
sheet_name <- if ("Normalized" %in% all_sheets) {
  "Normalized"
} else {
  all_sheets[1]
}
cat("Using sheet:", sheet_name, "\n\n")

# ── 3. Read data ──────────────────────────────────────────────────────────────
# Layout (CSF data.xlsx convention):
#   Row 1 : group labels  — first 3 cols = metadata, cols 4+ = "PD"/"Healthy"/etc.
#   Row 2 : column headers
#   Row 3+: protein rows  (log2-normalized intensities)

# Step A — read row 1 for group labels (no header so col positions are exact)
row1 <- read_excel(INPUT_FILE, sheet = sheet_name,
                   col_names = FALSE, n_max = 1)
n_cols_total <- ncol(row1)

cat(sprintf("Total columns in sheet: %d\n", n_cols_total))
cat("Row-1 unique values:", paste(unique(as.character(unlist(row1))), collapse = " | "), "\n\n")

# All values in row 1 (full row, not just cols 4+)
row1_all <- as.character(unlist(row1[1, ], use.names = FALSE))

# Normalise group labels regardless of position
LABEL_MAP <- c(
  "Control" = "Healthy", "HC" = "Healthy", "Normal" = "Healthy",
  "healthy" = "Healthy", "control" = "Healthy", "hc" = "Healthy",
  "PD" = "PD", "Parkinson" = "PD", "pd" = "PD", "parkinson" = "PD"
)
row1_norm <- dplyr::recode(row1_all, !!!LABEL_MAP)

# Identify ALL sample columns (any column whose row-1 label is Healthy or PD)
valid_idx  <- which(row1_norm %in% c("Healthy", "PD"))
groups_vec <- row1_norm[valid_idx]
groups     <- factor(groups_vec, levels = c("Healthy", "PD"))

cat(sprintf("Groups detected — Healthy: %d  PD: %d  (out of %d total cols)\n",
            sum(groups == "Healthy"), sum(groups == "PD"), n_cols_total))

if (length(valid_idx) == 0) {
  stop(paste0(
    "No group labels ('Healthy'/'Control'/'HC'/'PD') found in row 1.\n",
    "Row-1 contents: ", paste(row1_all, collapse = " | ")
  ))
}

# Step B — protein data (row 2 → header)
raw <- read_excel(INPUT_FILE, sheet = sheet_name, skip = 1)
cat(sprintf("Proteins loaded: %d  |  Columns after skip: %d\n", nrow(raw), ncol(raw)))

# Metadata: col 1 = Accession, col 2 = Gene Symbol (positional)
accessions   <- as.character(raw[[1]])
gene_symbols <- as.character(raw[[2]])

# Map valid_idx (1-based in row1) → column names in raw
# raw has same column order as row1 (read_excel adds header from row 2 = skip 1)
sample_cols <- names(raw)[valid_idx]
stopifnot(length(sample_cols) == length(groups_vec))

cat(sprintf("Sample columns used: %d  |  Group labels: %d\n\n",
            length(sample_cols), length(groups_vec)))

# ── 4. Expression matrix ──────────────────────────────────────────────────────
mat <- as.matrix(raw[, sample_cols])
storage.mode(mat) <- "numeric"
rownames(mat) <- seq_len(nrow(mat))

cat(sprintf("NA summary: %d proteins with ≥1 NA  |  %d fully complete\n",
            sum(rowSums(is.na(mat)) > 0),
            sum(rowSums(is.na(mat)) == 0)))

# ── NA filter: keep proteins observed in ≥50% of samples in EACH group ────────
min_obs_frac <- 0.50
hc_idx <- which(groups == "Healthy")
pd_idx <- which(groups == "PD")
min_hc <- ceiling(length(hc_idx) * min_obs_frac)
min_pd <- ceiling(length(pd_idx) * min_obs_frac)

keep <- (rowSums(!is.na(mat[, hc_idx, drop = FALSE])) >= min_hc) &
        (rowSums(!is.na(mat[, pd_idx, drop = FALSE])) >= min_pd)

cat(sprintf("After ≥50%% per-group filter: %d / %d proteins retained\n",
            sum(keep), nrow(mat)))

mat       <- mat[keep, ]
acc_keep  <- accessions[keep]
gene_keep <- gene_symbols[keep]

# Impute remaining NAs with the per-protein minimum observed value / 2
# (standard left-censored / MNAR approach for proteomics)
for (i in seq_len(nrow(mat))) {
  na_pos <- is.na(mat[i, ])
  if (any(na_pos)) {
    mat[i, na_pos] <- min(mat[i, !na_pos], na.rm = TRUE) / 2
  }
}

cat(sprintf("PON1 in set: %s | PON2 in set: %s\n\n",
            ifelse("PON1" %in% gene_keep, "YES", "NO"),
            ifelse("PON2" %in% gene_keep, "YES", "NO")))

# ── 5. limma DE analysis ──────────────────────────────────────────────────────
design   <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

cont_mat <- makeContrasts(PD_vs_Healthy = PD - Healthy, levels = design)
fit      <- lmFit(mat, design)
fit2     <- contrasts.fit(fit, cont_mat)
fit2     <- eBayes(fit2, trend = TRUE, robust = TRUE)

tt <- topTable(fit2, coef = "PD_vs_Healthy",
               number = Inf, sort.by = "none")
stopifnot(nrow(tt) == length(gene_keep))

results <- tibble(
  Accession = acc_keep,
  Gene      = gene_keep,
  log2FC    = tt$logFC,
  AveExpr   = tt$AveExpr,
  t_stat    = tt$t,
  pval      = tt$P.Value,
  adj_pval  = tt$adj.P.Val,
  B         = tt$B
) %>%
  mutate(
    sig       = adj_pval < 0.05 & abs(log2FC) > 0.58,
    direction = case_when(
      sig & log2FC > 0 ~ "UP",
      sig & log2FC < 0 ~ "DOWN",
      TRUE             ~ "NS"
    ),
    is_PON1 = Gene == "PON1",
    is_PON2 = Gene == "PON2",
    is_PON  = Gene %in% c("PON1", "PON2")
  )

cat(sprintf("limma results: UP=%d | DOWN=%d (adj.P<0.05, |log2FC|>0.58)\n",
            sum(results$direction == "UP"),
            sum(results$direction == "DOWN")))

for (g in c("PON1", "PON2")) {
  r <- filter(results, Gene == g)
  if (nrow(r) == 0) {
    cat(sprintf("  %s: not detected (filtered as NA)\n", g))
  } else {
    cat(sprintf("  %s | log2FC=%+.3f | p=%.4f | adj.p=%.4f | %s\n",
                g, r$log2FC, r$pval, r$adj_pval,
                ifelse(r$log2FC < 0, "DOWN in PD", "UP in PD")))
  }
}
cat("\n")

write_csv(results, file.path(OUTPUT_DIR, "limma_results_all.csv"))
cat("Saved: limma_results_all.csv\n\n")

# ── 6. Visual theme & helpers ─────────────────────────────────────────────────
COL_HEALTHY <- "#4575B4"
COL_PD      <- "#D73027"
COL_PON1    <- "#E69F00"
COL_PON2    <- "#FF4500"
COL_UP      <- "#CC4444"
COL_DOWN    <- "#4477AA"
COL_NS      <- "grey70"

BASE_THEME <- theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.title    = element_text(face = "bold"),
    legend.position = "none"
  )

fmt_pval <- function(p) {
  if (p < 0.001) "p < 0.001"
  else if (p < 0.01) "p < 0.01"
  else sprintf("p = %.3f", p)
}
fmt_adjp <- function(p) {
  if (p < 0.001) "adj.p < 0.001"
  else if (p < 0.01) "adj.p < 0.01"
  else sprintf("adj.p = %.3f", p)
}

n_healthy <- sum(groups == "Healthy")
n_pd      <- sum(groups == "PD")

# ── Plot 1 : PON1 + PON2 Boxplot (side-by-side facets) ───────────────────────
cat("[Plot 1] PON1 & PON2 boxplot\n")

pon_genes_present <- intersect(c("PON1", "PON2"), gene_keep)

if (length(pon_genes_present) == 0) {
  cat("  Neither PON1 nor PON2 detected — boxplot skipped.\n")
} else {
  # Build long-format expression table
  pon_long <- map_dfr(pon_genes_present, function(g) {
    idx <- which(gene_keep == g)[1]
    tibble(
      Gene       = g,
      Expression = as.numeric(mat[idx, ]),
      Group      = groups
    )
  })

  # Per-gene y positions for significance brackets
  pon_annot <- results %>%
    filter(Gene %in% pon_genes_present) %>%
    select(Gene, pval, adj_pval, log2FC) %>%
    group_by(Gene) %>%
    mutate(
      y_max   = max(pon_long$Expression[pon_long$Gene == Gene], na.rm = TRUE),
      rng     = diff(range(pon_long$Expression[pon_long$Gene == Gene], na.rm = TRUE)),
      y_seg   = y_max + rng * 0.08,
      y_lbl   = y_max + rng * 0.18,
      y_lim   = y_max + rng * 0.48,
      label   = paste0(sapply(pval,     fmt_pval), "\n",
                       sapply(adj_pval, fmt_adjp)),
      x_label = 1.5
    ) %>%
    ungroup()

  # Recompute per-gene stats outside dplyr to avoid row-wise issues
  pon_annot <- results %>%
    filter(Gene %in% pon_genes_present) %>%
    select(Gene, pval, adj_pval, log2FC) %>%
    rowwise() %>%
    mutate(
      gene_expr = list(pon_long$Expression[pon_long$Gene == Gene]),
      y_max     = max(unlist(gene_expr), na.rm = TRUE),
      rng       = diff(range(unlist(gene_expr), na.rm = TRUE)),
      y_seg     = y_max + rng * 0.08,
      y_lbl     = y_max + rng * 0.20,
      y_lim     = y_max + rng * 0.50,
      label     = paste0(fmt_pval(pval), "\n", fmt_adjp(adj_pval)),
      x_label   = 1.5
    ) %>%
    ungroup() %>%
    select(-gene_expr)

  p1 <- ggplot(pon_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.85,
                 color = "grey25", linewidth = 0.65) +
    geom_jitter(aes(color = Group), width = 0.12, size = 2.0, alpha = 0.70) +
    # significance bracket & label (per facet via blank + geom_text + segment)
    geom_blank(data = pon_annot, aes(x = 1, y = y_lim), inherit.aes = FALSE) +
    geom_segment(data = pon_annot,
                 aes(x = 1, xend = 2, y = y_seg, yend = y_seg),
                 inherit.aes = FALSE, linewidth = 0.8, color = "black") +
    geom_segment(data = pon_annot,
                 aes(x = 1, xend = 1,
                     y = y_max + (y_seg - y_max) * 0.25, yend = y_seg),
                 inherit.aes = FALSE, linewidth = 0.8, color = "black") +
    geom_segment(data = pon_annot,
                 aes(x = 2, xend = 2,
                     y = y_max + (y_seg - y_max) * 0.25, yend = y_seg),
                 inherit.aes = FALSE, linewidth = 0.8, color = "black") +
    geom_text(data = pon_annot,
              aes(x = x_label, y = y_lbl, label = label),
              inherit.aes = FALSE,
              size = 3.8, fontface = "bold", lineheight = 1.4) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = c(Healthy = COL_HEALTHY, PD = COL_PD),
      labels = c(Healthy = sprintf("Healthy (n=%d)", n_healthy),
                 PD      = sprintf("PD (n=%d)",      n_pd))
    ) +
    scale_color_manual(
      values = c(Healthy = COL_HEALTHY, PD = COL_PD)
    ) +
    scale_x_discrete(
      labels = c(Healthy = sprintf("Healthy\n(n=%d)", n_healthy),
                 PD      = sprintf("PD\n(n=%d)",      n_pd))
    ) +
    labs(
      title    = "PON1 & PON2 Expression in Substantia Nigra — Healthy vs PD",
      subtitle = "limma | log2 Normalized Intensity",
      x        = NULL,
      y        = expression(log[2] ~ "Normalized Intensity"),
      fill     = "Group"
    ) +
    BASE_THEME +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold", size = 13),
      strip.background = element_rect(fill = "grey92", color = NA)
    )

  ggsave(file.path(OUTPUT_DIR, "plot1_PON1_PON2_boxplot.pdf"), p1,
         width = 8, height = 7)
  ggsave(file.path(OUTPUT_DIR, "plot1_PON1_PON2_boxplot.png"), p1,
         width = 8, height = 7, dpi = 180)
  cat("  Saved: plot1_PON1_PON2_boxplot.pdf/.png\n")
}

# ── Plot 2 : Volcano plot — all proteins, PON1 & PON2 highlighted ─────────────
cat("[Plot 2] Volcano plot with PON1 & PON2 highlighted\n")

vol_dat <- results %>%
  mutate(
    log10p    = -log10(pmax(pval, 1e-300)),
    dot_color = case_when(
      is_PON1             ~ "PON1",
      is_PON2             ~ "PON2",
      direction == "UP"   ~ "UP",
      direction == "DOWN" ~ "DOWN",
      TRUE                ~ "NS"
    ),
    dot_size  = case_when(is_PON1 | is_PON2 ~ 5.5, TRUE ~ 1.8),
    dot_alpha = case_when(is_PON1 | is_PON2 ~ 1.0, TRUE ~ 0.50)
  ) %>%
  arrange(is_PON)   # draw PON proteins last (on top)

color_scale <- c(
  PON1 = COL_PON1,
  PON2 = COL_PON2,
  UP   = COL_UP,
  DOWN = COL_DOWN,
  NS   = COL_NS
)

# Build per-gene annotation labels
pon_vol <- filter(vol_dat, is_PON) %>%
  mutate(
    stat_label = paste0(
      Gene,
      "\nlog2FC = ", sprintf("%+.3f", log2FC), "\n",
      sapply(pval,     fmt_pval), "\n",
      sapply(adj_pval, fmt_adjp)
    )
  )

# Nudge directions so labels don't overlap each other
nudge_df <- pon_vol %>%
  mutate(
    nx = ifelse(log2FC < 0, -0.8, 0.8),
    ny = 1.5
  )

p2 <- ggplot(vol_dat, aes(x = log2FC, y = log10p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey45", linewidth = 0.55) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed",
             color = "grey45", linewidth = 0.55) +
  # Background proteins
  geom_point(data = filter(vol_dat, !is_PON),
             aes(color = dot_color, size = dot_size, alpha = dot_alpha)) +
  # PON1 & PON2 on top
  geom_point(data = filter(vol_dat, is_PON),
             aes(color = dot_color), size = 5.5, alpha = 1.0, shape = 18) +
  # Labels for PON1 & PON2
  geom_label_repel(
    data          = filter(nudge_df, Gene == "PON1"),
    aes(label     = stat_label),
    size          = 3.8, fontface = "bold",
    fill          = "#FFF9E6", color = COL_PON1,
    box.padding   = 1.0, point.padding = 0.6,
    segment.color = COL_PON1, segment.size  = 0.7,
    nudge_x = filter(nudge_df, Gene == "PON1")$nx,
    nudge_y = filter(nudge_df, Gene == "PON1")$ny,
    lineheight = 1.4, max.overlaps = Inf
  ) +
  geom_label_repel(
    data          = filter(nudge_df, Gene == "PON2"),
    aes(label     = stat_label),
    size          = 3.8, fontface = "bold",
    fill          = "#FFF3E0", color = COL_PON2,
    box.padding   = 1.0, point.padding = 0.6,
    segment.color = COL_PON2, segment.size  = 0.7,
    nudge_x = filter(nudge_df, Gene == "PON2")$nx,
    nudge_y = filter(nudge_df, Gene == "PON2")$ny,
    lineheight = 1.4, max.overlaps = Inf
  ) +
  scale_color_manual(
    values = color_scale,
    breaks = c("UP", "DOWN", "PON1", "PON2", "NS"),
    labels = c(
      UP   = sprintf("UP in PD (n=%d)",   sum(results$direction == "UP")),
      DOWN = sprintf("DOWN in PD (n=%d)", sum(results$direction == "DOWN")),
      PON1 = "PON1",
      PON2 = "PON2",
      NS   = "Not significant"
    ),
    name = NULL
  ) +
  scale_size_identity() +
  scale_alpha_identity() +
  annotate("text", x = -Inf, y = Inf,
           label = "DOWN in PD", hjust = -0.1, vjust = 1.8,
           color = "grey45", size = 3.5, fontface = "italic") +
  annotate("text", x =  Inf, y = Inf,
           label = "UP in PD",   hjust = 1.1, vjust = 1.8,
           color = "grey45", size = 3.5, fontface = "italic") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title    = "Volcano Plot — SN Proteomics (PD vs Healthy)",
    subtitle = sprintf("limma | %d proteins | PON1 & PON2 highlighted", nrow(results)),
    x        = expression(log[2] ~ "Fold Change (PD / Healthy)"),
    y        = expression(-log[10] ~ italic(P) * "-value")
  ) +
  BASE_THEME +
  theme(legend.position = "right", legend.text = element_text(size = 10))

ggsave(file.path(OUTPUT_DIR, "plot2_volcano_PON1_PON2.pdf"), p2,
       width = 9, height = 6.5)
ggsave(file.path(OUTPUT_DIR, "plot2_volcano_PON1_PON2.png"), p2,
       width = 9, height = 6.5, dpi = 180)
cat("  Saved: plot2_volcano_PON1_PON2.pdf/.png\n")

cat(sprintf("\nAll outputs saved to: %s/\n", OUTPUT_DIR))
cat(strrep("=", 55), "\n")
