# =============================================================================
# Dataset : PXD037684 (mmc1.xlsx — Supplementary Data)
# Study   : Mass Spectrometry-Based Proteomics Analysis of Human
#           Substantia Nigra From Parkinson's Disease Patients
# Journal : Molecular & Cellular Proteomics 22(1):100452, 2023
# Method  : Orbitrap + 11-plex TMT (3 batches) | PD n=15 vs HC n=15
#
# 정규화 방법 (논문과 동일):
#   1. 각 샘플값 ÷ 배치 내 MP(master pool)  → 배치 간 스케일 보정
#   2. log2 변환
#   3. ComBat                                → 잔여 배치 효과 제거
#   4. limma (~ 0 + group + batch)          → 차등발현 분석
# =============================================================================

# ── 0. 파일 경로 설정 ─────────────────────────────────────────────────────────
INPUT_FILE  <- "C:/Users/haaaa/Documents/mmc1.xlsx"
OUTPUT_DIR  <- "C:/Users/haaaa/Documents/PON_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. 패키지 로드 ────────────────────────────────────────────────────────────
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

cat("=== PXD037684 | PON1 & PON2 — Real Data Analysis ===\n")
cat("Input :", INPUT_FILE, "\n\n")

# ── 2. 데이터 읽기 ────────────────────────────────────────────────────────────
raw <- read_excel(INPUT_FILE, sheet = "ProteinTable")
cat(sprintf("Loaded: %d proteins x %d columns\n", nrow(raw), ncol(raw)))

# ── 3. Abundance 컬럼 추출 ────────────────────────────────────────────────────
# 컬럼명 예: "Abundances (Grouped): F1, 126_HC"
abund_cols <- names(raw)[grepl("^Abundances \\(Grouped\\)", names(raw))]

hc_cols <- abund_cols[grepl("_HC$", abund_cols)]   # 15개
pd_cols <- abund_cols[grepl("_PD$", abund_cols)]   # 15개
mp_cols <- abund_cols[grepl("_MP$", abund_cols)]   # 3개 (배치별 master pool)

cat(sprintf("Abundance columns: HC=%d, PD=%d, MP=%d\n",
            length(hc_cols), length(pd_cols), length(mp_cols)))

# 배치 레이블 추출 (F1/F2/F3)
get_batch <- function(col_names) {
  sub(".*: (F[123]),.*", "\\1", col_names)
}

# ── 4. 배치별 MP 정규화 → log2 변환 ──────────────────────────────────────────
# 논문 Methods: "intensity values were divided by the MP included in each set"
gene_col <- "Gene Symbol"
genes    <- raw[[gene_col]]

norm_list <- list()
for (b in c("F1", "F2", "F3")) {
  hc_b <- hc_cols[grepl(paste0(b, ","), hc_cols)]
  pd_b <- pd_cols[grepl(paste0(b, ","), pd_cols)]
  mp_b <- mp_cols[grepl(paste0(b, ","), mp_cols)]

  mp_vals <- raw[[mp_b]]                    # MP 값 (벡터, 길이=단백질수)
  mp_vals[mp_vals == 0] <- NA

  for (col in c(hc_b, pd_b)) {
    v <- raw[[col]] / mp_vals               # sample / MP
    v[v == 0] <- NA
    norm_list[[col]] <- v
  }
}

mat_norm <- do.call(cbind, norm_list)       # 단백질 x 30샘플
rownames(mat_norm) <- genes

# log2 변환
mat_log <- log2(mat_norm)

# 결측값이 하나라도 있는 단백질 제거 (논문과 동일 기준)
keep_rows <- rowSums(is.na(mat_log)) == 0
mat_log   <- mat_log[keep_rows, ]
cat(sprintf("After NA removal: %d proteins retained\n", nrow(mat_log)))

# ── 5. ComBat 배치 보정 ───────────────────────────────────────────────────────
# 샘플 메타데이터
sample_cols   <- c(hc_cols, pd_cols)
group_labels  <- ifelse(grepl("_HC$", sample_cols), "HC", "PD")
batch_labels  <- get_batch(sample_cols)

# ComBat: 배치를 공변량으로, 그룹을 보호
mod      <- model.matrix(~ group_labels)
mat_cb   <- ComBat(
  dat     = mat_log[, sample_cols],
  batch   = batch_labels,
  mod     = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

cat("ComBat batch correction done.\n")

# ── 6. limma 차등발현 분석 ────────────────────────────────────────────────────
groups <- factor(group_labels, levels = c("HC", "PD"))
batch  <- factor(batch_labels)
design <- model.matrix(~ 0 + groups + batch)
colnames(design) <- make.names(gsub("groups|batch", "", colnames(design)))

cont_mat <- makeContrasts(PD_vs_HC = PD - HC, levels = design)
fit  <- lmFit(mat_cb, design)
fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

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
    is_PON1 = Gene == "PON1",
    is_PON2 = Gene == "PON2",
    is_PON  = Gene %in% c("PON1", "PON2")
  )

cat(sprintf("DEP: UP %d | DOWN %d (adj.P<0.05, |log2FC|>0.58)\n\n",
            sum(results$direction == "UP"),
            sum(results$direction == "DOWN")))

# PON1/PON2 결과 출력
for (g in c("PON1", "PON2")) {
  r <- filter(results, Gene == g)
  if (nrow(r) == 0) {
    cat(sprintf("  %s : 데이터 없음 (결측값으로 제거됨)\n", g))
  } else {
    cat(sprintf("  %s | log2FC=%+.3f | p=%.3f | adj.p=%.3f | %s\n",
                g, r$log2FC, r$pval, r$adj_pval,
                ifelse(r$log2FC < 0, "DOWN in PD", "UP in PD")))
  }
}
cat("\n")

# ── 7. 시각화 설정 ────────────────────────────────────────────────────────────
COL_HC   <- "#4575B4"
COL_PD   <- "#D73027"
COL_PON1 <- "#E69F00"
COL_PON2 <- "#FF4500"

BASE_THEME <- theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.title    = element_text(face = "bold"),
    legend.position = "none"
  )

# p-value 포맷 함수
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

# ── Plot 1. PON2 Boxplot ──────────────────────────────────────────────────────
cat("[Plot 1] PON2 boxplot\n")

if ("PON2" %in% results$Gene) {
  pon2_mat <- mat_cb["PON2", sample_cols]
  pon2_dat <- tibble(
    Expression = as.numeric(pon2_mat),
    Group      = factor(group_labels, levels = c("HC", "PD")),
    Batch      = batch_labels
  )
  pon2_row <- filter(results, Gene == "PON2")
  y_max    <- max(pon2_dat$Expression, na.rm = TRUE)
  y_top    <- y_max + 1.2   # y축 상한: 2줄 annotation이 잘리지 않도록 여유 확보

  # p-value와 adj.p-value 두 줄로 표시
  pval_label <- fmt_pval(as.numeric(pon2_row$pval))
  adjp_label <- fmt_adjp(as.numeric(pon2_row$adj_pval))
  sig_label  <- paste0(pval_label, "\n", adjp_label)

  p1 <- ggplot(pon2_dat, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.85,
                 color = "grey25", linewidth = 0.65) +
    geom_jitter(aes(color = Group, shape = Batch),
                width = 0.12, size = 2.5, alpha = 0.80) +
    annotate("segment",
             x = 1, xend = 2, y = y_max + 0.20, yend = y_max + 0.20,
             linewidth = 0.8, color = "black") +
    annotate("segment",
             x = 1, xend = 1, y = y_max + 0.08, yend = y_max + 0.20,
             linewidth = 0.8, color = "black") +
    annotate("segment",
             x = 2, xend = 2, y = y_max + 0.08, yend = y_max + 0.20,
             linewidth = 0.8, color = "black") +
    annotate("text",
             x = 1.5, y = y_max + 0.45,
             label = sig_label,
             size = 3.8, fontface = "bold", lineheight = 1.4) +
    scale_fill_manual(values  = c(HC = COL_HC, PD = COL_PD)) +
    scale_color_manual(values = c(HC = COL_HC, PD = COL_PD)) +
    scale_shape_manual(values = c(F1 = 16, F2 = 17, F3 = 15),
                       name   = "TMT Batch") +
    scale_x_discrete(labels = c(HC = "HC\n(n=15)", PD = "PD\n(n=15)")) +
    scale_y_continuous(limits = c(NA, y_top)) +
    guides(shape = guide_legend(title = "TMT Batch")) +
    labs(
      title    = "PON2 Expression in Substantia Nigra",
      subtitle = "PXD037684 | Mol Cell Proteomics 2023 | Real Data",
      x        = NULL,
      y        = expression(log[2] ~ "Normalized TMT Intensity (sample/MP)")
    ) +
    BASE_THEME + theme(legend.position = "right")

  ggsave(file.path(OUTPUT_DIR, "plot1_PON2_boxplot.pdf"), p1,
         width = 5.5, height = 6.5)
  ggsave(file.path(OUTPUT_DIR, "plot1_PON2_boxplot.png"), p1,
         width = 5.5, height = 6.5, dpi = 180)
  cat("  Saved: plot1_PON2_boxplot.pdf/.png\n")
} else {
  cat("  PON2 skipped (no data after NA filter)\n")
}

# ── Plot 2. PON1 + PON2 Side-by-Side Boxplot ──────────────────────────────────
cat("[Plot 2] PON1 + PON2 boxplot\n")

pon_genes_present <- intersect(c("PON1", "PON2"), rownames(mat_cb))

if (length(pon_genes_present) > 0) {
  pon_long <- map_dfr(pon_genes_present, function(g) {
    tibble(
      Gene       = g,
      Expression = as.numeric(mat_cb[g, sample_cols]),
      Group      = factor(group_labels, levels = c("HC", "PD")),
      Batch      = batch_labels
    )
  })

  pon_ymax <- pon_long %>%
    group_by(Gene) %>%
    summarise(y = max(Expression, na.rm = TRUE) + 0.15, .groups = "drop")

  pon_pvals <- results %>%
    filter(Gene %in% pon_genes_present) %>%
    select(Gene, pval, adj_pval, log2FC) %>%
    left_join(pon_ymax, by = "Gene") %>%
    mutate(
      pval_label = sapply(pval,     fmt_pval),
      adjp_label = sapply(adj_pval, fmt_adjp),
      label      = paste0(pval_label, "\n", adjp_label),
      y_top      = y + 1.0   # 각 패널의 y 상한 (2줄 annotation 여유 포함)
    )

  p2 <- ggplot(pon_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.85,
                 color = "grey25", linewidth = 0.65) +
    geom_jitter(aes(color = Group), width = 0.12, size = 2.0, alpha = 0.75) +
    # 투명 더미 포인트로 각 패널 y축 상한 확장 (free_y 스케일 대응)
    geom_blank(data = pon_pvals, aes(x = 1, y = y_top),
               inherit.aes = FALSE) +
    geom_text(data = pon_pvals,
              aes(x = 1.5, y = y + 0.38, label = label),
              inherit.aes = FALSE, size = 3.5, fontface = "bold",
              lineheight = 1.4) +
    geom_segment(data = pon_pvals,
                 aes(x = 1, xend = 2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(HC = COL_HC, PD = COL_PD),
                      labels = c(HC = "HC (n=15)", PD = "PD (n=15)")) +
    scale_color_manual(values = c(HC = COL_HC, PD = COL_PD)) +
    scale_x_discrete(labels  = c(HC = "HC\n(n=15)", PD = "PD\n(n=15)")) +
    labs(
      title    = "PON1 & PON2 Expression in Substantia Nigra",
      subtitle = "PXD037684 | Mol Cell Proteomics 2023 | Real Data",
      x        = NULL,
      y        = expression(log[2] ~ "Normalized TMT Intensity (sample/MP)"),
      fill     = "Group"
    ) +
    BASE_THEME +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold", size = 13),
      strip.background = element_rect(fill = "grey92", color = NA)
    )

  ggsave(file.path(OUTPUT_DIR, "plot2_PON1_PON2_boxplot.pdf"), p2,
         width = 8, height = 7.0)
  ggsave(file.path(OUTPUT_DIR, "plot2_PON1_PON2_boxplot.png"), p2,
         width = 8, height = 7.0, dpi = 180)
  cat("  Saved: plot2_PON1_PON2_boxplot.pdf/.png\n")
}

# ── Plot 3. Volcano Plot (PON1 & PON2 only) ───────────────────────────────────
cat("[Plot 3] Volcano plot (PON1 & PON2 only)\n")

pon_dat <- results %>%
  filter(is_PON) %>%
  mutate(
    log10p     = -log10(pmax(pval, 1e-10)),
    pval_label = sapply(pval,     fmt_pval),
    adjp_label = sapply(adj_pval, fmt_adjp),
    stat_label = paste0(Gene,
                        "\nlog2FC = ", sprintf("%+.3f", log2FC),
                        "\n", pval_label,
                        "\n", adjp_label)
  )

if (nrow(pon_dat) == 0) {
  cat("  PON1/PON2 모두 결측값으로 제거됨 — volcano 생략\n")
} else {
  p3 <- ggplot(pon_dat, aes(x = log2FC, y = log10p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey45", linewidth = 0.55) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed",
               color = "grey45", linewidth = 0.55) +
    geom_point(data = filter(pon_dat, Gene == "PON1"),
               color = COL_PON1, size = 6, shape = 18) +
    geom_point(data = filter(pon_dat, Gene == "PON2"),
               color = COL_PON2, size = 6, shape = 18) +
    geom_label_repel(
      data          = filter(pon_dat, Gene == "PON1"),
      aes(label     = stat_label),
      size = 3.8, fontface = "bold",
      fill = "#FFF9E6", color = COL_PON1,
      box.padding = 1.0, point.padding = 0.6,
      segment.color = COL_PON1, segment.size = 0.6,
      nudge_x = -0.5, nudge_y = 1.0,
      lineheight = 1.4, max.overlaps = Inf
    ) +
    geom_label_repel(
      data          = filter(pon_dat, Gene == "PON2"),
      aes(label     = stat_label),
      size = 3.8, fontface = "bold",
      fill = "#FFF3E0", color = COL_PON2,
      box.padding = 1.0, point.padding = 0.6,
      segment.color = COL_PON2, segment.size = 0.6,
      nudge_x = 0.5, nudge_y = 1.0,
      lineheight = 1.4, max.overlaps = Inf
    ) +
    annotate("text", x = -Inf, y = Inf,
             label = "DOWN in PD", hjust = -0.1, vjust = 1.8,
             color = "grey50", size = 3.8, fontface = "bold") +
    annotate("text", x =  Inf, y = Inf,
             label = "UP in PD", hjust = 1.1, vjust = 1.8,
             color = "grey50", size = 3.8, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.25))) +
    labs(
      title    = "PON1 & PON2 — Substantia Nigra Proteomics (PD vs. HC)",
      subtitle = "PXD037684 | Real Data | 11-plex TMT | n=15/group",
      x = expression(log[2] ~ "Fold Change (PD / HC)"),
      y = expression(-log[10] ~ italic(P) - value)
    ) +
    BASE_THEME

  ggsave(file.path(OUTPUT_DIR, "plot3_volcano_PON1_PON2.pdf"), p3,
         width = 7, height = 6)
  ggsave(file.path(OUTPUT_DIR, "plot3_volcano_PON1_PON2.png"), p3,
         width = 7, height = 6, dpi = 180)
  cat("  Saved: plot3_volcano_PON1_PON2.pdf/.png\n")
}

# ── 8. 결과 CSV 저장 ──────────────────────────────────────────────────────────
write_csv(results, file.path(OUTPUT_DIR, "limma_results_all.csv"))
cat("\nAll results saved to:", OUTPUT_DIR, "\n")
cat(strrep("=", 55), "\n")
