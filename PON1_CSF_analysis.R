# =============================================================================
# Dataset : CSF data.xlsx (Normalized sheet)
# Study   : Parkinson's Disease CSF Proteomics
# Method  : log2-normalized CSF proteomics | PD n=40 vs Control n=40
#
# Analysis:
#   1. limma differential expression (PD vs Healthy/Control)
#   2. PON1 boxplot with p-value annotation
#   3. Volcano plot with PON1 highlighted
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma", ask = FALSE)
  for (pkg in c("readxl", "tidyverse", "ggrepel", "ggpubr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg)
  }
  library(limma)
  library(readxl)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
})

# ── 1. File paths ─────────────────────────────────────────────────────────────
INPUT_FILE <- "CSF data.xlsx"
OUTPUT_DIR <- "PON1_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== PD CSF Proteomics | PON1 Expression Analysis ===\n")
cat("Input:", INPUT_FILE, "\n\n")

# ── 2. Read data ──────────────────────────────────────────────────────────────
# Sheet layout:
#   Row 1 : group labels — cols 1-3 empty, cols 4-83 = PD / Control
#   Row 2 : headers     — Accession, Gene symbol, p-value, CSF01..CSF80
#   Row 3+ : protein data (log2-normalized values)
#
# Strategy: read the whole sheet without any headers to keep full position
# control, then extract group labels and expression values by numeric index.

raw_all <- read_excel(INPUT_FILE, sheet = "Normalized",
                      col_names = FALSE, col_types = "text")

n_rows <- nrow(raw_all)
n_cols <- ncol(raw_all)
cat(sprintf("Sheet dimensions: %d rows x %d columns\n", n_rows, n_cols))

# Row 1 = group labels; row 2 = header; rows 3+ = protein data
# Expression columns start at column 4 (0-based: index 4 in R = col D)
EXPR_START <- 4L
expr_col_idx <- EXPR_START:n_cols   # 4:83 → 80 columns

# Group labels from row 1, expression columns only
group_labels <- as.character(unlist(raw_all[1, expr_col_idx],
                                    use.names = FALSE))
group_labels[group_labels == "Control"] <- "Healthy"
groups <- factor(group_labels, levels = c("Healthy", "PD"))

cat(sprintf("Group labels extracted: %d  (Healthy=%d, PD=%d)\n",
            length(group_labels),
            sum(groups == "Healthy"), sum(groups == "PD")))

# Protein metadata from rows 3+, columns 1-3
accessions   <- as.character(unlist(raw_all[3:n_rows, 1], use.names = FALSE))
gene_symbols <- as.character(unlist(raw_all[3:n_rows, 2], use.names = FALSE))
sample_names <- as.character(unlist(raw_all[2, expr_col_idx], use.names = FALSE))

cat(sprintf("Proteins: %d | Sample columns: %d\n",
            length(accessions), length(sample_names)))

# Gene symbol lookup (Accession → Gene)
gene_map <- tibble(Accession = accessions, Gene = gene_symbols)

# ── 3. Build expression matrix ────────────────────────────────────────────────
# Extract numeric expression values; rows = proteins, cols = samples
mat_raw <- raw_all[3:n_rows, expr_col_idx]
mat <- matrix(
  as.numeric(unlist(mat_raw, use.names = FALSE)),
  nrow = length(accessions),
  ncol = length(sample_names),
  byrow = FALSE
)
rownames(mat) <- accessions      # unique Accession IDs as row identifiers
colnames(mat) <- sample_names

# Verify PON1 before NA filtering
pon1_acc <- gene_map$Accession[gene_map$Gene == "PON1"]
cat(sprintf("PON1 Accession: %s | NAs in row: %d\n",
            paste(pon1_acc, collapse = ","),
            sum(is.na(mat[pon1_acc, ]))))

# Remove rows with any NA
keep <- rowSums(is.na(mat)) == 0
mat  <- mat[keep, ]
gene_map <- gene_map[gene_map$Accession %in% rownames(mat), ]

cat(sprintf("After NA removal: %d proteins retained\n", nrow(mat)))
cat(sprintf("PON1 retained: %s\n\n",
            ifelse(any(pon1_acc %in% rownames(mat)), "YES", "NO")))

# ── 4. limma DE analysis ──────────────────────────────────────────────────────
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

cont_mat <- makeContrasts(PD_vs_Healthy = PD - Healthy, levels = design)
fit  <- lmFit(mat, design)
fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- topTable(fit2, coef = "PD_vs_Healthy", number = Inf,
                    sort.by = "none") %>%
  rownames_to_column("Accession") %>%
  as_tibble() %>%
  left_join(gene_map, by = "Accession") %>%       # attach Gene symbol
  rename(log2FC = logFC, pval = P.Value, adj_pval = adj.P.Val) %>%
  mutate(
    sig       = adj_pval < 0.05 & abs(log2FC) > 0.58,
    direction = case_when(
      sig & log2FC > 0 ~ "UP",
      sig & log2FC < 0 ~ "DOWN",
      TRUE             ~ "NS"
    ),
    is_PON1 = Gene == "PON1"
  )

cat(sprintf("\nlimma results: UP=%d | DOWN=%d (adj.P<0.05, |log2FC|>0.58)\n",
            sum(results$direction == "UP"),
            sum(results$direction == "DOWN")))

# PON1 result summary
pon1_row <- filter(results, Gene == "PON1")
if (nrow(pon1_row) == 0) {
  stop("PON1 not found after NA filtering — check input data.")
}
cat(sprintf("\nPON1 | log2FC=%+.3f | p=%.4f | adj.p=%.4f | %s\n\n",
            pon1_row$log2FC, pon1_row$pval, pon1_row$adj_pval,
            ifelse(pon1_row$log2FC < 0, "DOWN in PD", "UP in PD")))

# Save full limma results
write_csv(results, file.path(OUTPUT_DIR, "limma_results_all.csv"))
cat("Saved: limma_results_all.csv\n\n")

# ── 5. Color / theme settings ──────────────────────────────────────────────────
COL_HEALTHY <- "#4575B4"
COL_PD      <- "#D73027"
COL_PON1    <- "#E69F00"
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
  if (p < 0.001) "p < 0.001" else
    if (p < 0.01)  "p < 0.01"  else
      sprintf("p = %.3f", p)
}
fmt_adjp <- function(p) {
  if (p < 0.001) "adj.p < 0.001" else
    if (p < 0.01)  "adj.p < 0.01"  else
      sprintf("adj.p = %.3f", p)
}

# ── Plot 1. PON1 Boxplot ──────────────────────────────────────────────────────
cat("[Plot 1] PON1 boxplot (Healthy vs PD)\n")

pon1_expr <- as.numeric(mat[pon1_acc[1], ])
pon1_dat  <- tibble(
  Expression = pon1_expr,
  Group      = groups
)

y_max <- max(pon1_dat$Expression, na.rm = TRUE)
y_seg <- y_max + diff(range(pon1_dat$Expression, na.rm = TRUE)) * 0.08
y_lbl <- y_seg + diff(range(pon1_dat$Expression, na.rm = TRUE)) * 0.10
y_lim <- y_lbl + diff(range(pon1_dat$Expression, na.rm = TRUE)) * 0.22

pval_str <- fmt_pval(pon1_row$pval)
adjp_str <- fmt_adjp(pon1_row$adj_pval)
sig_label <- paste0(pval_str, "\n", adjp_str)

p1 <- ggplot(pon1_dat, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(
    width = 0.45, outlier.shape = NA, alpha = 0.85,
    color = "grey25", linewidth = 0.65
  ) +
  geom_jitter(
    aes(color = Group),
    width = 0.12, size = 2.0, alpha = 0.65
  ) +
  # significance bracket
  annotate("segment",
           x = 1, xend = 2,
           y = y_seg, yend = y_seg,
           linewidth = 0.8, color = "black") +
  annotate("segment",
           x = 1, xend = 1,
           y = y_max + diff(range(pon1_dat$Expression, na.rm = TRUE)) * 0.01,
           yend = y_seg,
           linewidth = 0.8, color = "black") +
  annotate("segment",
           x = 2, xend = 2,
           y = y_max + diff(range(pon1_dat$Expression, na.rm = TRUE)) * 0.01,
           yend = y_seg,
           linewidth = 0.8, color = "black") +
  annotate("text",
           x = 1.5, y = y_lbl,
           label = sig_label,
           size = 4.0, fontface = "bold", lineheight = 1.4) +
  scale_fill_manual(values = c(Healthy = COL_HEALTHY, PD = COL_PD)) +
  scale_color_manual(values = c(Healthy = COL_HEALTHY, PD = COL_PD)) +
  scale_x_discrete(
    labels = c(Healthy = sprintf("Healthy\n(n=%d)", sum(groups == "Healthy")),
               PD      = sprintf("PD\n(n=%d)",      sum(groups == "PD")))
  ) +
  coord_cartesian(ylim = c(NA, y_lim)) +
  labs(
    title    = "PON1 Expression in CSF — Healthy vs PD",
    subtitle = "limma | log2 Normalized Intensity",
    x        = NULL,
    y        = expression(log[2] ~ "Normalized Intensity")
  ) +
  BASE_THEME

ggsave(file.path(OUTPUT_DIR, "plot1_PON1_boxplot.pdf"), p1,
       width = 5, height = 6.5)
ggsave(file.path(OUTPUT_DIR, "plot1_PON1_boxplot.png"), p1,
       width = 5, height = 6.5, dpi = 180)
cat("  Saved: plot1_PON1_boxplot.pdf/.png\n")

# ── Plot 2. Volcano Plot (all proteins, PON1 highlighted) ─────────────────────
cat("[Plot 2] Volcano plot with PON1 highlighted\n")

vol_dat <- results %>%
  mutate(
    log10p    = -log10(pmax(pval, 1e-300)),
    dot_color = case_when(
      is_PON1             ~ "PON1",
      direction == "UP"   ~ "UP",
      direction == "DOWN" ~ "DOWN",
      TRUE                ~ "NS"
    ),
    dot_size  = if_else(is_PON1, 5, 1.8),
    dot_alpha = if_else(is_PON1, 1.0, 0.55)
  ) %>%
  arrange(is_PON1)  # draw PON1 on top

color_scale <- c(
  PON1 = COL_PON1,
  UP   = COL_UP,
  DOWN = COL_DOWN,
  NS   = COL_NS
)

pon1_vol <- filter(vol_dat, is_PON1)
pon1_label <- sprintf(
  "PON1\nlog2FC = %+.3f\n%s\n%s",
  pon1_vol$log2FC,
  fmt_pval(pon1_vol$pval),
  fmt_adjp(pon1_vol$adj_pval)
)

p2 <- ggplot(vol_dat, aes(x = log2FC, y = log10p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey45", linewidth = 0.55) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed",
             color = "grey45", linewidth = 0.55) +
  geom_point(
    aes(color = dot_color, size = dot_size, alpha = dot_alpha)
  ) +
  geom_label_repel(
    data          = pon1_vol,
    aes(label     = pon1_label),
    size          = 3.8,
    fontface      = "bold",
    fill          = "#FFF9E6",
    color         = COL_PON1,
    box.padding   = 1.2,
    point.padding = 0.7,
    segment.color = COL_PON1,
    segment.size  = 0.7,
    nudge_x       = -1.5,
    nudge_y       = 1.5,
    lineheight    = 1.4,
    max.overlaps  = Inf
  ) +
  scale_color_manual(
    values = color_scale,
    breaks = c("UP", "DOWN", "PON1", "NS"),
    labels = c(
      UP   = sprintf("UP in PD (n=%d)", sum(results$direction == "UP")),
      DOWN = sprintf("DOWN in PD (n=%d)", sum(results$direction == "DOWN")),
      PON1 = "PON1",
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
           label = "UP in PD", hjust = 1.1, vjust = 1.8,
           color = "grey45", size = 3.5, fontface = "italic") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title    = "Volcano Plot — CSF Proteomics (PD vs Healthy)",
    subtitle = sprintf("limma | %d proteins | PON1 highlighted",
                       nrow(results)),
    x = expression(log[2] ~ "Fold Change (PD / Healthy)"),
    y = expression(-log[10] ~ italic(P) * "-value")
  ) +
  BASE_THEME +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 10)
  )

ggsave(file.path(OUTPUT_DIR, "plot2_volcano_PON1.pdf"), p2,
       width = 8, height = 6.5)
ggsave(file.path(OUTPUT_DIR, "plot2_volcano_PON1.png"), p2,
       width = 8, height = 6.5, dpi = 180)
cat("  Saved: plot2_volcano_PON1.pdf/.png\n")

cat(sprintf("\nAll outputs saved to: %s/\n", OUTPUT_DIR))
cat(strrep("=", 55), "\n")
