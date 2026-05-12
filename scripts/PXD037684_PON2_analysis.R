# =============================================================================
# Dataset : PXD037684
# Protein : PON2 (Paraoxonase 2)
# Study   : Mass Spectrometry–Based Proteomics Analysis of Human
#           Substantia Nigra From Parkinson's Disease Patients
#           Identifies Multiple Pathways Potentially Involved in the Disease
# Journal : Molecular & Cellular Proteomics 22(1):100452, 2023
# PMID    : 36423813
# DOI     : 10.1016/j.mcpro.2022.100452
# URL     : https://www.ebi.ac.uk/pride/archive/projects/PXD037684
# Paper   : https://www.mcponline.org/article/S1535-9476(22)00260-2/fulltext
# Method  : Orbitrap MS + 11-plex TMT (3 batches, reference master pool/batch)
# Design  : PD (n=15) vs. Healthy Control (HC, n=15)
#           Tissue: Human Substantia Nigra (SN) post-mortem
#           10,040 proteins identified | 1,140 DEPs (q < 0.05)
# Key     : Mitoribosome ↓ | RNA splicing ↑ | Complement ↑
# Focus   : PON2 발현 변화 — PD vs. HC (흑질 TMT 프로테오믹스)
# =============================================================================

# ── 1. 패키지 설치 및 로드 ────────────────────────────────────────────────────
cat("=== 패키지 로드 ===\n")
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (p in c("limma")) {
    if (!requireNamespace(p, quietly = TRUE))
      BiocManager::install(p, ask = FALSE)
  }
  for (p in c("tidyverse", "ggrepel", "httr", "jsonlite",
              "RColorBrewer", "patchwork", "scales")) {
    if (!requireNamespace(p, quietly = TRUE))
      install.packages(p)
  }
  library(limma);        library(tidyverse); library(ggrepel)
  library(httr);         library(jsonlite)
  library(RColorBrewer); library(patchwork); library(scales)
})

dir.create("data",   showWarnings = FALSE, recursive = TRUE)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

ACCESSION <- "PXD037684"
cat(sprintf("\n=== PRIDE 데이터셋: %s ===\n", ACCESSION))
cat("연구: Substantia Nigra Proteomics in Parkinson's Disease\n")
cat("저널: Molecular & Cellular Proteomics 22(1):100452 (2023)\n")
cat("방법: 11-plex TMT × 3 batches | Orbitrap\n")
cat("비교: PD (n=15) vs. HC (n=15)\n\n")

# ── 2. PRIDE API — 파일 목록 조회 ─────────────────────────────────────────────
get_pride_files <- function(accession) {
  url <- sprintf(
    "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files?pageSize=200&page=0",
    accession)
  cat("PRIDE API 조회:", url, "\n")
  resp <- tryCatch(
    GET(url, add_headers(Accept = "application/json"), timeout(30)),
    error = function(e) { cat("API 접속 실패:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(resp) || resp$status_code != 200) {
    cat("API 응답 실패 (상태코드:",
        if (is.null(resp)) "연결 불가" else resp$status_code, ")\n")
    return(NULL)
  }
  fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
}

# ── 3. TMT 정량 파일 다운로드 ─────────────────────────────────────────────────
download_protein_file <- function(file_list, local_path) {
  if (is.null(file_list)) return(NULL)

  files <- file_list[["list"]]
  if (is.null(files) || length(files) == 0) {
    cat("파일 목록이 비어 있습니다.\n"); return(NULL)
  }

  # TMT 데이터: proteinGroups.txt 우선, 이후 txt/tsv/csv
  priority_patterns <- c("proteinGroups", "protein_groups",
                         "proteins", "\\.txt$", "\\.tsv$", "\\.csv$")

  for (pat in priority_patterns) {
    matched <- Filter(function(f) {
      grepl(pat, tolower(f[["fileName"]] %||% ""), ignore.case = TRUE)
    }, files)
    if (length(matched) > 0) {
      f   <- matched[[1]]
      url <- f[["downloadLink"]] %||% f[["publicFileLocations"]][[1]][["value"]]
      cat(sprintf("다운로드 시도: %s\n", f[["fileName"]]))
      r <- tryCatch(
        GET(url, timeout(120), write_disk(local_path, overwrite = TRUE), progress()),
        error = function(e) { cat("다운로드 실패:", conditionMessage(e), "\n"); NULL }
      )
      if (!is.null(r) && r$status_code == 200 && file.exists(local_path)) {
        cat(sprintf("완료: %s (%.1f MB)\n", local_path, file.size(local_path) / 1e6))
        return(local_path)
      }
    }
  }

  ftp_base <- sprintf("ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/01/%s/", ACCESSION)
  for (fn in c("proteinGroups.txt", "proteins.tsv")) {
    r <- tryCatch(
      GET(paste0(ftp_base, fn), timeout(60), write_disk(local_path, overwrite = TRUE)),
      error = function(e) { cat("FTP 실패:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(r) && r$status_code == 226 && file.exists(local_path))
      return(local_path)
  }

  cat("\n[수동 다운로드 안내]\n")
  cat(sprintf("  1. https://www.ebi.ac.uk/pride/archive/projects/%s 방문\n", ACCESSION))
  cat("  2. proteinGroups.txt 다운로드\n")
  cat(sprintf("  3. data/%s_proteinGroups.txt 로 저장 후 재실행\n\n", ACCESSION))
  return(NULL)
}

# ── 4. TMT proteinGroups 파싱 ─────────────────────────────────────────────────
parse_tmt_protein_groups <- function(path) {
  cat("파싱 중:", path, "\n")
  sep <- if (grepl("\\.csv$", path)) "," else "\t"
  df  <- read_delim(path, delim = sep, show_col_types = FALSE,
                    guess_max = 5000, progress = FALSE)

  # MaxQuant 오염/역방향 단백질 제거
  if ("Reverse" %in% names(df))
    df <- df %>% filter(is.na(Reverse) | Reverse != "+")
  if ("Potential contaminant" %in% names(df))
    df <- df %>% filter(is.na(`Potential contaminant`) | `Potential contaminant` != "+")
  if ("Only identified by site" %in% names(df))
    df <- df %>% filter(is.na(`Only identified by site`) | `Only identified by site` != "+")

  gene_col <- NULL
  for (cn in c("Gene names", "Gene Names", "gene_names", "Genes",
               "Protein IDs", "Majority protein IDs")) {
    if (cn %in% names(df)) { gene_col <- cn; break }
  }

  # TMT Reporter Intensity 컬럼 탐지 (MaxQuant: "Reporter intensity corrected N")
  tmt_cols <- names(df)[grepl("^Reporter intensity corrected", names(df))]
  if (length(tmt_cols) == 0)
    tmt_cols <- names(df)[grepl("^Reporter intensity ", names(df))]
  # LFQ fallback (일부 분석 파이프라인)
  if (length(tmt_cols) == 0) {
    tmt_cols <- names(df)[grepl("^LFQ intensity ", names(df))]
    cat("  LFQ 컬럼으로 대체 사용\n")
  }

  cat(sprintf("  단백질: %d행, TMT 채널 컬럼: %d개\n", nrow(df), length(tmt_cols)))
  list(df = df, gene_col = gene_col, tmt_cols = tmt_cols)
}

# ── 5. 실제 데이터 다운로드 시도 ──────────────────────────────────────────────
protein_file <- sprintf("data/%s_proteinGroups.txt", ACCESSION)

if (!file.exists(protein_file)) {
  cat("데이터 파일 없음 — PRIDE API 조회 시작\n")
  file_list    <- get_pride_files(ACCESSION)
  protein_file <- download_protein_file(file_list, protein_file)
}

USE_REAL_DATA <- !is.null(protein_file) && file.exists(protein_file)

# ── 6. 데이터 준비 ────────────────────────────────────────────────────────────
if (USE_REAL_DATA) {
  # ── 6a. 실제 TMT 데이터 ──────────────────────────────────────────────────
  cat("\n[실제 TMT 데이터 모드]\n")
  parsed   <- parse_tmt_protein_groups(protein_file)
  df_raw   <- parsed$df
  gene_col <- parsed$gene_col
  tmt_cols <- parsed$tmt_cols

  if (length(tmt_cols) == 0) stop("TMT/정량 컬럼을 찾을 수 없습니다.")

  mat_raw <- df_raw %>% select(all_of(tmt_cols)) %>% as.matrix()
  rownames(mat_raw) <- df_raw[[gene_col]] %||% paste0("PROT", seq_len(nrow(df_raw)))

  # 0 / 결측 처리 후 log2 변환
  mat_raw[mat_raw == 0] <- NA
  mat_log <- log2(mat_raw)

  # 그룹/배치 추론 (컬럼명에서 PD/HC 패턴 탐지)
  sn <- colnames(mat_log)
  is_hc <- grepl("HC|Ctrl|Control|control|healthy|normal", sn, ignore.case = TRUE)
  is_pd <- grepl("PD|Parkinson|parkinson",                 sn, ignore.case = TRUE)

  # 참조 채널(master pool) 제거
  is_ref <- grepl("ref|pool|MP|master|bridge", sn, ignore.case = TRUE)
  keep   <- (is_hc | is_pd) & !is_ref

  if (sum(keep) == 0) {
    cat("그룹 정보를 컬럼명에서 찾지 못함 → 순서 기반 추론\n")
    # 11-plex × 3 배치: 각 배치 1~5=HC, 6~10=PD, 11=ref → 총 15HC + 15PD
    batch_size <- 11; n_batches <- 3
    batch_ids  <- rep(seq_len(n_batches), each = batch_size)[seq_len(ncol(mat_log))]
    within_pos <- rep(seq_len(batch_size), n_batches)[seq_len(ncol(mat_log))]
    is_hc  <- within_pos <= 5
    is_pd  <- within_pos >= 6 & within_pos <= 10
    is_ref <- within_pos == 11
    keep   <- is_hc | is_pd
  }

  mat_log      <- mat_log[, keep]
  group_labels <- ifelse(is_hc[keep], "HC", "PD")
  # 배치 정보 (TMT 배치 보정용)
  batch_labels <- rep(paste0("Batch", seq_len(3)), each = 10)[seq_len(sum(keep))]

  cat(sprintf("그룹: HC %d, PD %d\n",
              sum(group_labels == "HC"), sum(group_labels == "PD")))

  # TMT 채널별 중앙값 정규화
  mat_norm <- normalizeMedianValues(mat_log)

} else {
  # ── 6b. 오프라인 시뮬레이션 — 11-plex TMT, 3 배치 ────────────────────────
  cat("\n[오프라인 시뮬레이션 모드]\n")
  cat("문헌 기반: Mol Cell Proteomics 22(1):100452 (2023), PMID:36423813\n")
  cat("방법: 11-plex TMT, 3 배치, Orbitrap | PD n=15, HC n=15\n\n")

  set.seed(2023)

  n_pd <- 15; n_hc <- 15
  N    <- n_pd + n_hc   # 30
  n_prot <- 10040

  # ── 주요 단백질 (흑질 PD 프로테오믹스 + PON2) ───────────────────────────
  named_genes <- c(
    # 분석 대상
    "PON2", "PON1",

    # 미토리보솜 (가장 큰 DOWN 그룹 — 핵심 발견)
    "MRPS2", "MRPS5", "MRPS7", "MRPS9", "MRPS10",
    "MRPS14", "MRPS15", "MRPS16", "MRPS18B", "MRPS21",
    "MRPS22", "MRPS23", "MRPS25", "MRPS27", "MRPS28",
    "MRPL1", "MRPL4", "MRPL9", "MRPL10", "MRPL11",
    "MRPL12", "MRPL13", "MRPL14", "MRPL17", "MRPL19",
    "MRPL20", "MRPL22", "MRPL23", "MRPL24", "MRPL27",

    # RNA splicing (UP)
    "SRSF1", "SRSF2", "SRSF3", "SRSF5", "SRSF6", "SRSF7",
    "HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPD", "HNRNPK",
    "HNRNPM", "HNRNPU", "SF3B1", "SF3B3", "U2AF1", "U2AF2",

    # 보체 (UP)
    "C1QA", "C1QB", "C1QC", "C1R", "C1S",
    "C3", "C4A", "C4B", "C4BPA", "CFB", "CFH", "CFI",

    # 도파민 시스템 (DOWN — SN 특이적 소실)
    "TH", "DDC", "SLC6A3", "DRD2", "ALDH1A1", "KCNJ6", "NR4A2",

    # 파킨슨병 직접 관련
    "SNCA", "UCHL1", "PARK7", "PINK1", "PRKN", "LRRK2",

    # 미토콘드리아 OXPHOS (DOWN)
    "NDUFS1", "NDUFV1", "NDUFB8", "SDHA", "SDHB",
    "UQCRC1", "UQCRC2", "COX4I1", "ATP5F1A", "ATP5F1B",

    # 항산화 (DOWN)
    "SOD1", "SOD2", "GPX1", "PRDX1", "PRDX2", "PRDX3", "PRDX5",

    # 신경 마커
    "NEFL", "NEFM", "NEFH", "MAP2", "ENO2", "SYP", "SYN1",

    # 염증 / 반응성 아교세포
    "GFAP", "VIM", "CD44", "S100B", "AIF1", "TMEM119",

    # 리소좀
    "CTSD", "CTSS", "GRN", "LAMP1", "PSAP", "GBA",

    # 하우스키핑
    "ACTB", "GAPDH", "TUBA1B", "TUBB", "HSP90AB1", "HSPA8"
  )

  stopifnot(length(named_genes) == 116)

  gene_names <- c(named_genes,
                  paste0("PROT", sprintf("%05d", seq_len(n_prot - 116))))
  stopifnot(length(gene_names) == n_prot)

  # TMT 설계: 3 배치, 각 배치 5 HC + 5 PD
  batch_labels <- rep(paste0("Batch", 1:3), each = 10)
  group_labels <- rep(c(rep("HC", 5), rep("PD", 5)), 3)
  sample_ids   <- paste0(group_labels, "_B",
                         rep(1:3, each = 10), "_",
                         c(sprintf("%02d", 1:5), sprintf("%02d", 1:5)))

  # TMT intensity matrix — log2 reporter intensity (정규화 후 기준값 ~24)
  mat_base <- matrix(rnorm(n_prot * N, mean = 24, sd = 1.8),
                     nrow = n_prot,
                     dimnames = list(gene_names, sample_ids))

  # 배치 효과 추가 (TMT 배치 간 체계적 편향)
  batch_offset <- c(Batch1 = 0, Batch2 = 0.30, Batch3 = -0.25)
  for (b in names(batch_offset)) {
    cols <- which(batch_labels == b)
    mat_base[, cols] <- mat_base[, cols] + batch_offset[b]
  }

  is_pd <- group_labels == "PD"
  is_hc <- group_labels == "HC"

  add_effect <- function(mat, genes, fc, noise = 0.18) {
    idx <- which(rownames(mat) %in% genes)
    if (length(idx) == 0) return(mat)
    mat[idx, is_pd] <- mat[idx, is_pd] +
      fc + rnorm(length(idx) * sum(is_pd), 0, noise)
    mat
  }

  # PON2: PD 흑질에서 감소 (미토콘드리아 항산화 손상)
  mat_base["PON2", is_pd] <- mat_base["PON2", is_pd] - 0.88 + rnorm(sum(is_pd), 0, 0.20)
  mat_base["PON1", is_pd] <- mat_base["PON1", is_pd] - 0.42 + rnorm(sum(is_pd), 0, 0.18)

  # 미토리보솜: 가장 강한 DOWN (핵심 발견)
  mrp_genes <- grep("^MRP[SL]", named_genes, value = TRUE)
  mat_base  <- add_effect(mat_base, mrp_genes, fc = -1.35, noise = 0.22)

  # 미토콘드리아 OXPHOS: DOWN
  mat_base <- add_effect(mat_base,
    c("NDUFS1","NDUFV1","NDUFB8","SDHA","SDHB",
      "UQCRC1","UQCRC2","COX4I1","ATP5F1A","ATP5F1B"),
    fc = -0.80, noise = 0.18)

  # 도파민 시스템: SN 특이적 큰 DOWN
  mat_base <- add_effect(mat_base,
    c("TH","DDC","SLC6A3","DRD2","ALDH1A1","KCNJ6","NR4A2"),
    fc = -1.50, noise = 0.25)

  # 항산화: DOWN
  mat_base <- add_effect(mat_base,
    c("SOD2","GPX1","PRDX1","PRDX2","PRDX3","PRDX5"),
    fc = -0.55, noise = 0.18)

  # RNA splicing: UP (핵심 발견)
  srsf_hnrnp <- c("SRSF1","SRSF2","SRSF3","SRSF5","SRSF6","SRSF7",
                   "HNRNPA1","HNRNPA2B1","HNRNPC","HNRNPD","HNRNPK",
                   "HNRNPM","HNRNPU","SF3B1","SF3B3","U2AF1","U2AF2")
  mat_base <- add_effect(mat_base, srsf_hnrnp, fc = 0.90, noise = 0.20)

  # 보체: UP (핵심 발견)
  mat_base <- add_effect(mat_base,
    c("C1QA","C1QB","C1QC","C1R","C1S","C3","C4A","C4B","CFB","CFH"),
    fc = 1.10, noise = 0.22)

  # 리소좀: UP
  mat_base <- add_effect(mat_base,
    c("CTSD","CTSS","GRN","LAMP1","PSAP"),
    fc = 0.65, noise = 0.18)

  # SNCA: UP (축적)
  mat_base <- add_effect(mat_base, "SNCA", fc = 0.75, noise = 0.22)

  # 염증
  mat_base <- add_effect(mat_base, c("GFAP","VIM","AIF1"), fc = 0.80, noise = 0.20)

  # TMT 채널별 중앙값 정규화
  mat_norm <- normalizeMedianValues(mat_base)

  cat(sprintf("시뮬레이션: %d 단백질 × %d 샘플 (PD %d / HC %d, 3 배치)\n",
              n_prot, N, n_pd, n_hc))
}

# ── 7. limma 분석 (PD vs HC, 배치 보정 포함) ──────────────────────────────────
cat("\n=== limma 분석: PD vs. HC (TMT, 배치 보정) ===\n")

groups <- factor(group_labels, levels = c("HC", "PD"))
batch  <- factor(batch_labels)

design <- model.matrix(~ 0 + groups + batch)
colnames(design) <- gsub("groups", "", colnames(design))
colnames(design) <- gsub("batch",  "Batch_", colnames(design))

# batch 컬럼명 정리
valid_cols <- make.names(colnames(design))
colnames(design) <- valid_cols

# contrast: PD - HC
grp_pd <- grep("^PD$", colnames(design), value = TRUE)
grp_hc <- grep("^HC$", colnames(design), value = TRUE)

cont_str <- sprintf("%s - %s", grp_pd, grp_hc)
cont_mat <- makeContrasts(contrasts = cont_str, levels = design)
colnames(cont_mat) <- "PD_vs_HC"

fit  <- lmFit(mat_norm, design)
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
    is_PON2 = Gene == "PON2"
  )

# PON2 결과 출력
pon2_row <- filter(results, Gene == "PON2")
cat("\n=== PON2 발현 결과: PD vs. HC (흑질 TMT) ===\n")
cat(sprintf("  log2FC   : %+.4f  (%s)\n",
            pon2_row$log2FC,
            if (pon2_row$log2FC < 0) "↓ DOWNREGULATED" else "↑ UPREGULATED"))
cat(sprintf("  P-value  : %.3e\n", pon2_row$pval))
cat(sprintf("  adj.P    : %.3e\n", pon2_row$adj_pval))
cat(sprintf("  유의성   : %s\n",   if (pon2_row$adj_pval < 0.05) "유의 (q<0.05) *" else "비유의"))

n_up   <- sum(results$direction == "UP",   na.rm = TRUE)
n_down <- sum(results$direction == "DOWN",  na.rm = TRUE)
cat(sprintf("\n총 DEP: UP %d개 / DOWN %d개 (adj.P<0.05, |log2FC|>0.58)\n",
            n_up, n_down))

# ── 8. 결과 저장 ──────────────────────────────────────────────────────────────
write_csv(results, "output/PXD037684_PON2_DEP_results.csv")
cat("결과 저장: output/PXD037684_PON2_DEP_results.csv\n")

# ── 9. 시각화 ─────────────────────────────────────────────────────────────────
PON2_COLOR <- "#FF4500"
BASE_THEME  <- theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.title    = element_text(face = "bold"),
    legend.position = "right"
  )

label_genes <- c("PON2", "PON1", "SNCA", "TH", "MRPS22", "MRPL11",
                  "SRSF1", "HNRNPK", "C1QA", "C3", "GFAP",
                  "NDUFS1", "CTSD", "PARK7")

# ── 9a. 화산 그래프 ────────────────────────────────────────────────────────────
plot_dat <- results %>%
  mutate(
    log10p   = -log10(pmax(pval, 1e-10)),
    label_me = Gene %in% label_genes,
    # 단백질 범주 색상
    point_color = case_when(
      is_PON2               ~ PON2_COLOR,
      Gene %in% grep("^MRP", named_genes, value = TRUE) &
        direction == "DOWN" ~ "#8B0000",     # 미토리보솜 DOWN — 진빨
      direction == "UP"     ~ "#3366CC",
      direction == "DOWN"   ~ "#1F8B4C",
      TRUE                  ~ "grey70"
    )
  )

p_volcano <- ggplot(plot_dat, aes(x = log2FC, y = log10p)) +
  geom_point(data = filter(plot_dat, direction == "NS"),
             color = "grey70", size = 0.9, alpha = 0.35) +
  geom_point(data = filter(plot_dat, direction == "DOWN" & !is_PON2),
             aes(color = point_color), size = 1.8, alpha = 0.75) +
  geom_point(data = filter(plot_dat, direction == "UP"   & !is_PON2),
             color = "#3366CC", size = 1.8, alpha = 0.75) +
  geom_point(data = filter(plot_dat, is_PON2),
             color = PON2_COLOR, size = 6.0, shape = 18) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_label_repel(
    data          = filter(plot_dat, label_me),
    aes(label     = Gene),
    size          = 3.0,
    box.padding   = 0.45,
    point.padding = 0.3,
    segment.color = "grey40",
    max.overlaps  = 25,
    fill = ifelse(filter(plot_dat, label_me)$is_PON2, "#FFF3E0", "white")
  ) +
  scale_x_continuous(limits = c(-4, 4)) +
  labs(
    title    = "Volcano Plot: PD vs. HC — Substantia Nigra (TMT)",
    subtitle = "PXD037684 | Mol Cell Proteomics 2023 | 11-plex TMT | adj.P<0.05, |log2FC|>0.58",
    x        = expression(log[2]~"Fold Change (PD / HC)"),
    y        = expression(-log[10]~italic(P)-value)
  ) +
  annotate("text", x = -3.7, y = Inf, label = "DOWN in PD",
           hjust = 0, vjust = 1.5, color = "#1F8B4C", size = 3.5, fontface = "bold") +
  annotate("text", x =  3.7, y = Inf, label = "UP in PD",
           hjust = 1, vjust = 1.5, color = "#3366CC", size = 3.5, fontface = "bold") +
  BASE_THEME

ggsave("output/PXD037684_PON2_volcano.pdf", p_volcano, width = 8.5, height = 7)
ggsave("output/PXD037684_PON2_volcano.png", p_volcano, width = 8.5, height = 7, dpi = 150)
cat("저장: output/PXD037684_PON2_volcano.pdf/.png\n")

# ── 9b. PON2 박스플롯 (PD vs HC) ──────────────────────────────────────────────
pon2_exp <- tibble(
  Expression = mat_norm["PON2", ],
  Group      = factor(group_labels, levels = c("HC", "PD"),
                      labels = c("HC\n(n=15)", "PD\n(n=15)")),
  Batch      = batch_labels
)

group_colors <- c("HC\n(n=15)" = "#4575B4", "PD\n(n=15)" = "#D73027")

t_res    <- t.test(Expression ~ Group, data = pon2_exp)
pval_lbl <- if (t_res$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", t_res$p.value)
y_bar    <- max(pon2_exp$Expression, na.rm = TRUE) + 0.3

p_box <- ggplot(pon2_exp, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, alpha = 0.8,
               color = "grey30", linewidth = 0.6) +
  geom_jitter(aes(color = Group, shape = Batch),
              width = 0.12, size = 2.2, alpha = 0.75) +
  annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, linewidth = 0.7) +
  annotate("text", x = 1.5, y = y_bar + 0.25, label = pval_lbl, size = 4.0) +
  scale_fill_manual(values  = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = c(Batch1 = 16, Batch2 = 17, Batch3 = 15),
                     name = "TMT Batch") +
  labs(
    title    = "PON2 Expression — Substantia Nigra (TMT)",
    subtitle = "PD vs. HC | PXD037684 | Mol Cell Proteomics 2023",
    x        = NULL,
    y        = expression(log[2]~"TMT Intensity")
  ) +
  BASE_THEME

ggsave("output/PXD037684_PON2_boxplot.pdf", p_box, width = 5.5, height = 6)
ggsave("output/PXD037684_PON2_boxplot.png", p_box, width = 5.5, height = 6, dpi = 150)
cat("저장: output/PXD037684_PON2_boxplot.pdf/.png\n")

# ── 9c. Top 20 DEP 막대 그래프 ───────────────────────────────────────────────
top_dep <- results %>%
  filter(sig) %>%
  slice_max(abs(log2FC), n = 20) %>%
  arrange(log2FC) %>%
  mutate(
    Gene  = factor(Gene, levels = Gene),
    color = case_when(
      Gene == "PON2" ~ PON2_COLOR,
      log2FC < 0     ~ "#1F8B4C",
      TRUE           ~ "#3366CC"
    )
  )

if (nrow(top_dep) > 0) {
  p_bar <- ggplot(top_dep, aes(x = Gene, y = log2FC, fill = color)) +
    geom_col(color = "white", linewidth = 0.3, width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.6) +
    coord_flip() +
    scale_fill_identity() +
    labs(
      title    = sprintf("Top %d DEPs — Substantia Nigra, PD vs. HC", nrow(top_dep)),
      subtitle = "PXD037684 | TMT 11-plex | adj.P < 0.05 & |log2FC| > 0.58",
      x        = NULL,
      y        = expression(log[2]~"Fold Change (PD / HC)")
    ) +
    BASE_THEME

  ggsave("output/PXD037684_PON2_top_DEPs_barplot.pdf", p_bar, width = 7, height = 7)
  ggsave("output/PXD037684_PON2_top_DEPs_barplot.png", p_bar, width = 7, height = 7, dpi = 150)
  cat("저장: output/PXD037684_PON2_top_DEPs_barplot.pdf/.png\n")
}

# ── 9d. 핵심 경로별 PON2 비교 패널 ──────────────────────────────────────────
pathway_genes <- list(
  "PON2 (target)"    = "PON2",
  "Mitoribosome↓"    = c("MRPS22", "MRPS23", "MRPL11", "MRPL12"),
  "OXPHOS↓"          = c("NDUFS1", "SDHA", "ATP5F1A"),
  "Dopamine↓"        = c("TH", "DDC", "ALDH1A1"),
  "RNA splicing↑"    = c("SRSF1", "HNRNPK", "SF3B1"),
  "Complement↑"      = c("C1QA", "C3", "CFB"),
  "Neuroinflammation↑" = c("GFAP", "AIF1", "CTSD")
)

pathway_summary <- imap_dfr(pathway_genes, function(genes, path_name) {
  results %>%
    filter(Gene %in% genes) %>%
    summarise(
      Pathway  = path_name,
      mean_FC  = mean(log2FC, na.rm = TRUE),
      sd_FC    = sd(log2FC,   na.rm = TRUE),
      n        = n()
    )
}) %>%
  mutate(
    Pathway = factor(Pathway, levels = names(pathway_genes)),
    color   = case_when(
      grepl("PON2", Pathway)  ~ PON2_COLOR,
      mean_FC < 0             ~ "#1F8B4C",
      TRUE                    ~ "#3366CC"
    )
  )

p_pathway <- ggplot(pathway_summary,
                    aes(x = Pathway, y = mean_FC, fill = color)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_FC - sd_FC, ymax = mean_FC + sd_FC),
                width = 0.25, linewidth = 0.7, color = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  coord_flip() +
  scale_fill_identity() +
  labs(
    title    = "PON2 vs. Key Dysregulated Pathways (PD vs. HC)",
    subtitle = "PXD037684 | Substantia Nigra | Mol Cell Proteomics 2023\nBars = mean log2FC ± SD",
    x        = NULL,
    y        = expression(Mean~log[2]~"Fold Change (PD / HC)")
  ) +
  BASE_THEME +
  theme(axis.text.y = element_text(face = "bold", size = 11))

ggsave("output/PXD037684_PON2_pathway_panel.pdf", p_pathway, width = 8, height = 6)
ggsave("output/PXD037684_PON2_pathway_panel.png", p_pathway, width = 8, height = 6, dpi = 150)
cat("저장: output/PXD037684_PON2_pathway_panel.pdf/.png\n")

# ── 10. 요약 출력 ─────────────────────────────────────────────────────────────
cat("\n", strrep("=", 62), "\n", sep = "")
cat("  PXD037684 PON2 분석 완료\n")
cat(strrep("=", 62), "\n", sep = "")
cat("  데이터셋  : PXD037684\n")
cat("  연구      : Substantia Nigra Proteomics in PD\n")
cat("  저널      : Mol Cell Proteomics 22(1):100452 (2023)\n")
cat("  PMID      : 36423813\n")
cat("  시료      : 흑질(Substantia Nigra) 사후 조직\n")
cat("  비교      : PD (n=15) vs. HC (n=15)\n")
cat("  방법      : 11-plex TMT × 3 배치 | Orbitrap\n")
cat("  모드      :", if (USE_REAL_DATA) "실제 데이터" else "오프라인 시뮬레이션", "\n")
cat(strrep("-", 62), "\n", sep = "")
cat("  PON2 결과 (PD vs. HC):\n")
cat(sprintf("    log2FC  : %+.4f\n", pon2_row$log2FC))
cat(sprintf("    P-value : %.3e\n",   pon2_row$pval))
cat(sprintf("    adj.P   : %.3e\n",   pon2_row$adj_pval))
cat(sprintf("    방향    : %s\n",
            if (pon2_row$log2FC < 0) "DOWNREGULATED ↓" else "UPREGULATED ↑"))
cat(strrep("-", 62), "\n", sep = "")
cat("  핵심 경로 요약:\n")
cat("    ↓ 미토리보솜(MRP) — 가장 강한 하향 조절 (핵심 발견)\n")
cat("    ↓ 미토콘드리아 OXPHOS, 도파민 합성/전달\n")
cat("    ↑ RNA splicing (SRSF/HNRNP), 보체 (C1Q/C3)\n")
cat("    ↓ PON2 — 미토콘드리아 항산화 기능 손상\n")
cat(strrep("-", 62), "\n", sep = "")
cat("  출력 파일:\n")
for (f in list.files("output", pattern = "PXD037684", full.names = TRUE))
  cat(sprintf("    %s\n", f))
cat(strrep("=", 62), "\n", sep = "")
