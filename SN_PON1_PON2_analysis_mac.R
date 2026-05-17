# =============================================================================
# Dataset : SN protein data.xlsx  (pon2 branch)
# Study   : Parkinson's Disease — Substantia Nigra Proteomics
# Format  : Thermo PD export | 11-plex TMT, 3 batches (F1/F2/F3)
#           Column names: "Abundances (Grouped): Fx, xxx_HC/PD/MP"
#
# Analysis:
#   1. Per-batch MP normalisation → log2
#   2. ComBat batch correction
#   3. limma DE (PD vs HC)
#   4. PON1 + PON2 side-by-side boxplot with p-value annotation
#   5. Full volcano plot with PON1 & PON2 highlighted
#
# macOS note: if the file is not found in the working directory, a
#             native Finder file-picker opens automatically.
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (pkg in c("limma", "sva")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE)
  }
  for (pkg in c("readxl", "tidyverse", "ggrepel", "ggpubr", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg)
  }
  library(limma);   library(sva)
  library(readxl);  library(tidyverse)
  library(ggrepel); library(ggpubr); library(patchwork)
})

# ── 1. File path — auto file-picker on macOS if not found ────────────────────
INPUT_FILE <- "SN protein data.xlsx"

if (!file.exists(INPUT_FILE)) {
  message("'SN protein data.xlsx' not found in working directory.")
  message("Opening file picker — please select the xlsx file.")
  INPUT_FILE <- file.choose()        # native macOS Finder dialog
  # Change working directory to the file location so output lands beside it
  setwd(dirname(INPUT_FILE))
}

OUTPUT_DIR <- file.path(dirname(INPUT_FILE), "SN_PON_output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== PD SN Proteomics | PON1 & PON2 Expression Analysis ===\n")
cat("Input :", INPUT_FILE, "\n")
cat("Output:", OUTPUT_DIR, "\n\n")

# ── 2. Sheet selection ────────────────────────────────────────────────────────
all_sheets <- excel_sheets(INPUT_FILE)
cat("Available sheets:", paste(all_sheets, collapse = ", "), "\n")
sheet_name <- all_sheets[1]
cat("Using sheet:", sheet_name, "\n\n")

# ── 3. Read data ──────────────────────────────────────────────────────────────
raw <- read_excel(INPUT_FILE, sheet = sheet_name)
cat(sprintf("Loaded: %d proteins x %d columns\n", nrow(raw), ncol(raw)))

# ── 4. Abundance column extraction ───────────────────────────────────────────
abund_cols <- names(raw)[grepl("^Abundances \\(Grouped\\)", names(raw))]
cat(sprintf("Abundance columns found: %d\n", length(abund_cols)))

hc_cols <- abund_cols[grepl("_HC$", abund_cols)]
pd_cols <- abund_cols[grepl("_PD$", abund_cols)]
mp_cols <- abund_cols[grepl("_MP$", abund_cols)]

cat(sprintf("  HC: %d  |  PD: %d  |  MP (master pool): %d\n\n",
            length(hc_cols), length(pd_cols), length(mp_cols)))

if (length(hc_cols) == 0 || length(pd_cols) == 0)
  stop("No _HC or _PD columns found. Check column name format.")

get_batch    <- function(cols) sub(".*: (F[0-9]+),.*", "\\1", cols)
batches_all  <- unique(get_batch(c(hc_cols, pd_cols)))
cat("Batches detected:", paste(sort(batches_all), collapse = ", "), "\n\n")

# ── 5. Gene symbol column ─────────────────────────────────────────────────────
gene_col_candidates <- c("Gene Symbol", "Gene.Symbol", "Gene", "gene_symbol",
                          "Genes", "gene", "GENE")
gene_col <- intersect(gene_col_candidates, names(raw))[1]
if (is.na(gene_col)) {
  cat("Available columns:\n"); print(names(raw)[1:min(30, ncol(raw))])
  stop("Gene symbol column not found. Edit gene_col_candidates above.")
}
cat("Gene column:", gene_col, "\n\n")

genes_raw <- as.character(raw[[gene_col]])
# Replace NA / empty entries before make.unique so matrix rownames are never NA
bad_gene          <- is.na(genes_raw) | trimws(genes_raw) == "" | genes_raw == "NA"
genes_raw[bad_gene] <- paste0("Unknown_", seq_along(genes_raw)[bad_gene])
genes             <- make.unique(genes_raw)

# ── 6. Per-batch MP normalisation → log2 ─────────────────────────────────────
norm_list <- list()

for (b in sort(batches_all)) {
  hc_b <- hc_cols[get_batch(hc_cols) == b]
  pd_b <- pd_cols[get_batch(pd_cols) == b]
  mp_b <- mp_cols[get_batch(mp_cols) == b]

  if (length(mp_b) == 0) {
    warning(sprintf("No MP column for batch %s; skipping MP normalisation.", b))
    for (col in c(hc_b, pd_b)) {
      v <- as.numeric(raw[[col]]); v[v <= 0] <- NA
      norm_list[[col]] <- log2(v)
    }
  } else {
    mp_vals <- as.numeric(raw[[mp_b[1]]]); mp_vals[mp_vals <= 0] <- NA
    for (col in c(hc_b, pd_b)) {
      v <- as.numeric(raw[[col]]) / mp_vals; v[v <= 0] <- NA
      norm_list[[col]] <- log2(v)
    }
  }
}

sample_cols  <- c(hc_cols, pd_cols)
group_labels <- c(rep("HC", length(hc_cols)), rep("PD", length(pd_cols)))
batch_labels <- get_batch(sample_cols)

mat_log <- do.call(cbind, norm_list[sample_cols])
rownames(mat_log) <- genes

cat(sprintf("Samples: HC=%d  PD=%d\n", sum(group_labels == "HC"),
            sum(group_labels == "PD")))
cat(sprintf("Proteins before NA filter: %d\n", nrow(mat_log)))

# ── 7. NA filter (≥50 % valid per group) ─────────────────────────────────────
hc_idx <- which(group_labels == "HC")
pd_idx <- which(group_labels == "PD")
min_hc <- ceiling(length(hc_idx) * 0.50)
min_pd <- ceiling(length(pd_idx) * 0.50)

keep_rows <- (rowSums(!is.na(mat_log[, hc_idx, drop = FALSE])) >= min_hc) &
             (rowSums(!is.na(mat_log[, pd_idx, drop = FALSE])) >= min_pd)
mat_log   <- mat_log[keep_rows, ]

cat(sprintf("After ≥50%% per-group filter: %d proteins retained\n", nrow(mat_log)))
cat(sprintf("PON1 retained: %s | PON2 retained: %s\n\n",
            ifelse(any(grepl("^PON1", rownames(mat_log))), "YES", "NO"),
            ifelse(any(grepl("^PON2", rownames(mat_log))), "YES", "NO")))

# ── 8. ComBat batch correction ────────────────────────────────────────────────
mat_imp <- mat_log
for (i in seq_len(nrow(mat_imp))) {
  na_pos <- is.na(mat_imp[i, ])
  if (any(na_pos))
    mat_imp[i, na_pos] <- min(mat_imp[i, !na_pos], na.rm = TRUE) / 2
}

mod    <- model.matrix(~ group_labels)
mat_cb <- ComBat(dat = mat_imp, batch = batch_labels, mod = mod,
                 par.prior = TRUE, prior.plots = FALSE)
cat("ComBat batch correction done.\n\n")

# ── 9. limma DE analysis ──────────────────────────────────────────────────────
groups <- factor(group_labels, levels = c("HC", "PD"))
batch  <- factor(batch_labels)
design <- model.matrix(~ 0 + groups + batch)
colnames(design) <- make.names(sub("groups|batch", "", colnames(design)))

cont_mat <- makeContrasts(PD_vs_HC = PD - HC, levels = design)
fit      <- lmFit(mat_cb, design)
fit2     <- contrasts.fit(fit, cont_mat)
fit2     <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- topTable(fit2, coef = "PD_vs_HC", number = Inf, sort.by = "none") %>%
  rownames_to_column("Gene") %>%
  as_tibble() %>%
  rename(log2FC = logFC, pval = P.Value, adj_pval = adj.P.Val) %>%
  mutate(
    sig       = adj_pval < 0.05 & abs(log2FC) > 0.58,
    direction = case_when(
      sig & log2FC > 0 ~ "UP",
      sig & log2FC < 0 ~ "DOWN",
      TRUE             ~ "NS"
    ),
    is_PON1 = grepl("^PON1(\\.\\d+)?$", Gene),
    is_PON2 = grepl("^PON2(\\.\\d+)?$", Gene),
    is_PON  = is_PON1 | is_PON2
  )

cat(sprintf("limma: UP=%d | DOWN=%d (adj.P<0.05, |log2FC|>0.58)\n\n",
            sum(results$direction == "UP"),
            sum(results$direction == "DOWN")))

for (g in c("PON1", "PON2")) {
  r <- filter(results, grepl(paste0("^", g, "(\\.\\d+)?$"), Gene)) %>%
    slice_min(pval, n = 1, with_ties = FALSE)
  if (nrow(r) == 0) {
    cat(sprintf("  %s: not detected (filtered out)\n", g))
  } else {
    cat(sprintf("  %s | log2FC=%+.3f | p=%.4f | adj.p=%.4f | %s\n",
                g, r$log2FC, r$pval, r$adj_pval,
                ifelse(r$log2FC < 0, "DOWN in PD", "UP in PD")))
  }
}
cat("\n")

write_csv(results, file.path(OUTPUT_DIR, "limma_results_all.csv"))
cat("Saved: limma_results_all.csv\n\n")

# ── 10. Visual theme & helpers ────────────────────────────────────────────────
COL_HC   <- "#4575B4"
COL_PD   <- "#D73027"
COL_PON1 <- "#E69F00"
COL_PON2 <- "#FF4500"
COL_UP   <- "#CC4444"
COL_DOWN <- "#4477AA"
COL_NS   <- "grey70"

BASE_THEME <- theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    axis.title      = element_text(face = "bold"),
    legend.position = "none"
  )

fmt_pval <- function(p) {
  if      (p < 0.001) "p < 0.001"
  else if (p < 0.01)  "p < 0.01"
  else                sprintf("p = %.3f", p)
}
fmt_adjp <- function(p) {
  if      (p < 0.001) "adj.p < 0.001"
  else if (p < 0.01)  "adj.p < 0.01"
  else                sprintf("adj.p = %.3f", p)
}

n_hc <- sum(group_labels == "HC")
n_pd <- sum(group_labels == "PD")

# ── Plot 1 : PON1 + PON2 side-by-side boxplot ─────────────────────────────────
cat("[Plot 1] PON1 & PON2 boxplot\n")

find_pon_rowname <- function(mat, display) {
  hits <- grep(paste0("^", display, "(\\.\\d+)?$"), rownames(mat), value = TRUE)
  if (length(hits) == 0) NULL else hits[1]
}
pon_map <- Filter(Negate(is.null), setNames(
  lapply(c("PON1", "PON2"), find_pon_rowname, mat = mat_cb),
  c("PON1", "PON2")
))

if (length(pon_map) == 0) {
  cat("  Neither PON1 nor PON2 detected — boxplot skipped.\n")
} else {
  pon_long <- map_dfr(names(pon_map), function(display) {
    rn <- pon_map[[display]]
    tibble(
      Gene       = display,
      Expression = as.numeric(mat_cb[rn, sample_cols]),
      Group      = factor(group_labels, levels = c("HC", "PD")),
      Batch      = batch_labels
    )
  })

  pon_annot <- map_dfr(names(pon_map), function(display) {
    rn  <- pon_map[[display]]
    r   <- filter(results, Gene == rn) %>%
           slice_min(pval, n = 1, with_ties = FALSE)
    if (nrow(r) == 0) return(tibble())
    dat   <- filter(pon_long, Gene == display)
    y_max <- max(dat$Expression, na.rm = TRUE)
    rng   <- diff(range(dat$Expression, na.rm = TRUE))
    tibble(
      Gene     = display,
      pval     = r$pval[1],     adj_pval = r$adj_pval[1],
      log2FC   = r$log2FC[1],   y_max    = y_max,
      y_seg    = y_max + rng * 0.08,
      y_lbl    = y_max + rng * 0.20,
      y_lim    = y_max + rng * 0.52,
      label    = paste0(fmt_pval(r$pval[1]), "\n", fmt_adjp(r$adj_pval[1])),
      x_label  = 1.5
    )
  })

  p1 <- ggplot(pon_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.85,
                 color = "grey25", linewidth = 0.65) +
    geom_jitter(aes(color = Group, shape = Batch),
                width = 0.12, size = 2.5, alpha = 0.80) +
    geom_blank(data = pon_annot, aes(x = 1, y = y_lim), inherit.aes = FALSE) +
    geom_segment(data = pon_annot,
                 aes(x = 1, xend = 2, y = y_seg, yend = y_seg),
                 inherit.aes = FALSE, linewidth = 0.8, color = "black") +
    geom_segment(data = pon_annot,
                 aes(x = 1, xend = 1, y = y_max + 0.02, yend = y_seg),
                 inherit.aes = FALSE, linewidth = 0.8, color = "black") +
    geom_segment(data = pon_annot,
                 aes(x = 2, xend = 2, y = y_max + 0.02, yend = y_seg),
                 inherit.aes = FALSE, linewidth = 0.8, color = "black") +
    geom_text(data = pon_annot, aes(x = x_label, y = y_lbl, label = label),
              inherit.aes = FALSE, size = 3.8, fontface = "bold",
              lineheight = 1.4) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = c(HC = COL_HC, PD = COL_PD),
      labels = c(HC = sprintf("HC (n=%d)", n_hc),
                 PD = sprintf("PD (n=%d)", n_pd))
    ) +
    scale_color_manual(values = c(HC = COL_HC, PD = COL_PD)) +
    scale_shape_manual(values = c(F1 = 16, F2 = 17, F3 = 15),
                       name   = "TMT Batch") +
    scale_x_discrete(
      labels = c(HC = sprintf("HC\n(n=%d)", n_hc),
                 PD = sprintf("PD\n(n=%d)", n_pd))
    ) +
    labs(
      title    = "PON1 & PON2 Expression in Substantia Nigra — HC vs PD",
      subtitle = "limma | log2(sample/MP) | ComBat batch correction",
      x        = NULL,
      y        = expression(log[2] ~ "(Normalised TMT Intensity)"),
      fill     = "Group"
    ) +
    BASE_THEME +
    theme(
      legend.position  = "right",
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
    dot_size  = if_else(is_PON1 | is_PON2, 5.5, 1.8),
    dot_alpha = if_else(is_PON1 | is_PON2, 1.0, 0.50)
  ) %>%
  arrange(is_PON)

color_scale <- c(
  PON1 = COL_PON1, PON2 = COL_PON2,
  UP   = COL_UP,   DOWN = COL_DOWN, NS = COL_NS
)

pon_vol <- filter(vol_dat, is_PON) %>%
  mutate(display = case_when(is_PON1 ~ "PON1", is_PON2 ~ "PON2")) %>%
  group_by(display) %>%
  slice_min(pval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stat_label = paste0(
      display,
      "\nlog2FC = ", sprintf("%+.3f", log2FC), "\n",
      map_chr(pval,     fmt_pval), "\n",
      map_chr(adj_pval, fmt_adjp)
    ),
    nx = if_else(log2FC < 0, -0.8, 0.8),
    ny = 1.5
  )

pon1_vol <- filter(pon_vol, display == "PON1")
pon2_vol <- filter(pon_vol, display == "PON2")

p2 <- ggplot(vol_dat, aes(x = log2FC, y = log10p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey45", linewidth = 0.55) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed",
             color = "grey45", linewidth = 0.55) +
  geom_point(data = filter(vol_dat, !is_PON),
             aes(color = dot_color, size = dot_size, alpha = dot_alpha)) +
  geom_point(data = filter(vol_dat, is_PON),
             aes(color = dot_color), size = 5.5, alpha = 1.0, shape = 18) +
  { if (nrow(pon1_vol) > 0) geom_label_repel(
      data = pon1_vol, aes(label = stat_label),
      size = 3.8, fontface = "bold",
      fill = "#FFF9E6", color = COL_PON1,
      box.padding = 1.0, point.padding = 0.6,
      segment.color = COL_PON1, segment.size = 0.7,
      nudge_x = pon1_vol$nx, nudge_y = pon1_vol$ny,
      lineheight = 1.4, max.overlaps = Inf
  ) else NULL } +
  { if (nrow(pon2_vol) > 0) geom_label_repel(
      data = pon2_vol, aes(label = stat_label),
      size = 3.8, fontface = "bold",
      fill = "#FFF3E0", color = COL_PON2,
      box.padding = 1.0, point.padding = 0.6,
      segment.color = COL_PON2, segment.size = 0.7,
      nudge_x = pon2_vol$nx, nudge_y = pon2_vol$ny,
      lineheight = 1.4, max.overlaps = Inf
  ) else NULL } +
  scale_color_manual(
    values = color_scale,
    breaks = c("UP", "DOWN", "PON1", "PON2", "NS"),
    labels = c(
      UP   = sprintf("UP in PD (n=%d)",   sum(results$direction == "UP")),
      DOWN = sprintf("DOWN in PD (n=%d)", sum(results$direction == "DOWN")),
      PON1 = "PON1", PON2 = "PON2",
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
    title    = "Volcano Plot — SN Proteomics (PD vs HC)",
    subtitle = sprintf("limma | %d proteins | PON1 & PON2 highlighted",
                       nrow(results)),
    x        = expression(log[2] ~ "Fold Change (PD / HC)"),
    y        = expression(-log[10] ~ italic(P) * "-value")
  ) +
  BASE_THEME +
  theme(legend.position = "right", legend.text = element_text(size = 10))

ggsave(file.path(OUTPUT_DIR, "plot2_volcano_PON1_PON2.pdf"), p2,
       width = 9, height = 6.5)
ggsave(file.path(OUTPUT_DIR, "plot2_volcano_PON1_PON2.png"), p2,
       width = 9, height = 6.5, dpi = 180)
cat("  Saved: plot2_volcano_PON1_PON2.pdf/.png\n")

cat(sprintf("\nAll outputs saved to:\n  %s\n", OUTPUT_DIR))
cat(strrep("=", 55), "\n")
