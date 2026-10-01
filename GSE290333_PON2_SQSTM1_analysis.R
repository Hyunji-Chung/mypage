# =============================================================================
# Dataset : GEO GSE290333 (Clariom S Human, PBMC, 23 de novo PD vs 16 CT)
# Paper   : Navarrete et al., Neurotherapeutics 22 (2025) e00762
#           "Whole transcriptome analysis of PBMCs from de novo and
#            drug-naive Parkinson's disease patients"
# Target  : PON2, SQSTM1
#
# Why     : The paper's limma model adjusted for sex only.
#           (1) PD (mean 70.6 y) vs CT (57.7 y) -> strong age confounding
#           (2) CIBERSORTx: PD had more neutrophils, fewer naive CD4 T cells
#           -> Is the small SQSTM1 increase (and PON2 decrease) a disease
#              effect, an age effect, or a cell-composition effect?
#
# Paper result (Supplementary mmc2, model ~ group + sex):
#   SQSTM1 TC0500009707.hg.1  logFC  0.383  P 3.5e-4  adj.P 0.110
#   SQSTM1 TC0500009706.hg.1  logFC  0.148  P 5.5e-3  adj.P 0.275
#   PON2   TC0700011835.hg.1  logFC -0.143  P 0.092   adj.P 0.556
#
# Analysis:
#   1. Download GSE290333 (series matrix; optional raw CEL + RMA)
#   2. Sample metadata (GEO characteristics, else paper Table 1/2)
#   3. Reproduce paper model (~ group + sex) and check against mmc2
#   4. Cell composition covariates (CIBERSORTx output if given,
#      else marker-gene scores)
#   5. Sequential limma models: + age, + cells, + age + cells, ...
#   6. Sensitivity: age-overlap subset, age outlier removed,
#      age slope within CT only / PD only
#   7. Plots: forest plot, age scatter, raw vs adjusted boxplot,
#      composition scores, correlation heatmap
# =============================================================================

# ── 0. Settings ───────────────────────────────────────────────────────────────
GSE_ID           <- "GSE290333"
DATA_DIR         <- "GSE290333_data"
OUTPUT_DIR       <- "GSE290333_PON2_SQSTM1_output"

# Supplementary Multimedia component 2 (paper DE table, 18,091 probes).
# Used for (a) reproducing the paper's probe filter, (b) probe -> gene
# annotation fallback, (c) checking the reproduction. Set NULL if absent.
PAPER_DE_FILE    <- "mmc2_5.xlsx"
FILTER_TO_PAPER_PROBES <- TRUE

# CIBERSORTx result file (CSV/TSV/XLSX; first column = sample name,
# e.g. "Mixture"). NULL -> marker-gene based cell scores are used instead.
CIBERSORT_FILE   <- NULL

# Optional manual metadata (CSV) — columns: geo_accession, group (CT/PD),
# sex (Female/Male), age [, site, rin]. Overrides everything else.
SAMPLE_META_FILE <- NULL

# TRUE -> re-normalize from raw CEL files (oligo::rma) instead of the
# series matrix. Needs pd.clariom.s.human.
USE_RAW_CEL      <- FALSE

TARGET_GENES  <- c("PON2", "SQSTM1")
TARGET_PROBES <- c("TC0500009707.hg.1",   # SQSTM1 (paper top hit)
                   "TC0500009706.hg.1",   # SQSTM1
                   "TC0700011835.hg.1")   # PON2

# ── 1. Packages ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (pkg in c("limma", "GEOquery", "Biobase")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  for (pkg in c("readxl", "tidyverse", "patchwork")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg)
  }
  library(limma)
  library(GEOquery)
  library(Biobase)
  library(readxl)
  library(tidyverse)
  library(patchwork)
})

options(timeout = max(600, getOption("timeout")))
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== GSE290333 PBMC | PON2 & SQSTM1 | age + cell-composition adjustment ===\n\n")

# ── 2. Paper metadata (Table 1 & 2) ──────────────────────────────────────────
# Fallback when GEO characteristics do not carry age/sex.
paper_meta <- tribble(
  ~sample_id, ~age, ~sex,     ~hy,   ~site,   ~rin,
  "CT1",  66, "Female", NA,    "12OUH", 7.3,
  "CT2",  60, "Female", NA,    "12OUH", 8.6,
  "CT3",  61, "Female", NA,    "GUHA",  7.4,
  "CT4",  58, "Male",   NA,    "GUHA",  8.4,
  "CT5",  66, "Female", NA,    "12OUH", 8.7,
  "CT6",  69, "Male",   NA,    "12OUH", 8.7,
  "CT7",  68, "Male",   NA,    "12OUH", 8.7,
  "CT8",  49, "Male",   NA,    "12OUH", 7.8,
  "CT9",  47, "Female", NA,    "12OUH", 8.8,
  "CT10", 49, "Female", NA,    "12OUH", 6.4,
  "CT11", 62, "Female", NA,    "12OUH", 6.5,
  "CT12", 50, "Female", NA,    "12OUH", 6.3,
  "CT13", 55, "Male",   NA,    "12OUH", 6.3,
  "CT14", 45, "Male",   NA,    "12OUH", 6.4,
  "CT15", 52, "Male",   NA,    "12OUH", 6.3,
  "CT16", 66, "Male",   NA,    "12OUH", 7.5,
  "PD1",  71, "Male",   "I",   "12OUH", 7.6,
  "PD2",  72, "Female", "II",  "12OUH", 7.7,
  "PD3",  62, "Female", "I",   "12OUH", 7.4,
  "PD4",  65, "Male",   "I",   "12OUH", 6.7,
  "PD5",  82, "Male",   "II",  "12OUH", 8.2,
  "PD6",  75, "Male",   "I",   "12OUH", 7.8,
  "PD7", 106, "Male",   "II",  "12OUH", 8.2,
  "PD8",  81, "Female", "II",  "12OUH", 8.4,
  "PD9",  74, "Male",   "I",   "12OUH", 8.1,
  "PD10", 71, "Male",   "I",   "GUHA",  9.7,
  "PD11", 69, "Male",   "II",  "GUHA",  9.6,
  "PD12", 70, "Female", "I",   "GUHA",  9.7,
  "PD13", 77, "Male",   "II",  "12OUH", 9.3,
  "PD14", 63, "Female", "I",   "12OUH", 5.1,
  "PD15", 71, "Female", "III", "12OUH", 7.3,
  "PD16", 84, "Female", "II",  "12OUH", 8.1,
  "PD17", 82, "Female", "II",  "12OUH", 8.9,
  "PD18", 78, "Female", "I",   "12OUH", 6.5,
  "PD19", 42, "Female", "I",   "GUHA",  9.4,
  "PD20", 55, "Female", "I",   "GUHA",  9.1,
  "PD21", 56, "Male",   "II",  "GUHA",  9.1,
  "PD22", 47, "Male",   "I",   "12OUH", 8.0,
  "PD23", 71, "Male",   "II",  "12OUH", 6.4
)

# ── 3. Expression data ───────────────────────────────────────────────────────
cat("[1] Downloading", GSE_ID, "series matrix ...\n")
gse   <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = TRUE, destdir = DATA_DIR)
eset  <- gse[[1]]
pheno <- pData(eset)
fdat  <- fData(eset)
expr  <- exprs(eset)

if (USE_RAW_CEL || nrow(expr) == 0 || all(is.na(expr))) {
  cat("    -> Using raw CEL files + RMA\n")
  for (pkg in c("oligo", "pd.clariom.s.human")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  getGEOSuppFiles(GSE_ID, baseDir = DATA_DIR, filter_regex = "RAW\\.tar")
  tar_file <- file.path(DATA_DIR, GSE_ID, paste0(GSE_ID, "_RAW.tar"))
  cel_dir  <- file.path(DATA_DIR, GSE_ID, "CEL")
  untar(tar_file, exdir = cel_dir)
  cels <- list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$",
                     ignore.case = TRUE, full.names = TRUE)
  raw  <- oligo::read.celfiles(cels)
  expr <- exprs(oligo::rma(raw, target = "core"))
  colnames(expr) <- str_extract(colnames(expr), "GSM[0-9]+")
  expr <- expr[, intersect(colnames(expr), rownames(pheno)), drop = FALSE]
}

if (max(expr, na.rm = TRUE) > 100) {
  cat("    -> Values look non-log; applying log2(x + 1)\n")
  expr <- log2(expr + 1)
}
expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]
cat("    Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n\n")

# ── 4. Probe annotation ───────────────────────────────────────────────────────
paper_de <- NULL
if (!is.null(PAPER_DE_FILE) && file.exists(PAPER_DE_FILE)) {
  paper_de <- read_excel(PAPER_DE_FILE) %>%
    rename(probe = Row.names) %>%
    select(probe, logFC, P.Value, adj.P.Val, GENE_NAME)
  cat("[2] Paper DE table loaded:", nrow(paper_de), "probes\n")
} else {
  cat("[2] Paper DE table not found (", PAPER_DE_FILE, ") — skipping checks\n")
}

symbol <- setNames(rep(NA_character_, nrow(expr)), rownames(expr))
if (nrow(fdat) > 0) {
  fd <- fdat[match(rownames(expr), rownames(fdat)), , drop = FALSE]
  direct <- intersect(c("GENE_SYMBOL", "Gene Symbol", "gene_symbol", "SYMBOL"),
                      colnames(fd))
  if (length(direct) > 0) {
    symbol[] <- as.character(fd[[direct[1]]])
  } else {
    # Affymetrix "gene_assignment" style: "NM_xxx // SYMBOL // desc // ..."
    ga_col <- colnames(fd)[vapply(fd, function(x)
      any(grepl(" // ", head(as.character(x), 2000))), logical(1))]
    if (length(ga_col) > 0) {
      symbol[] <- vapply(strsplit(as.character(fd[[ga_col[1]]]), " // "),
                         function(s) if (length(s) >= 2) trimws(s[2]) else NA_character_,
                         character(1))
    }
  }
}
symbol[symbol %in% c("", "---")] <- NA
if (!is.null(paper_de)) {
  miss <- is.na(symbol) & names(symbol) %in% paper_de$probe
  symbol[miss] <- paper_de$GENE_NAME[match(names(symbol)[miss], paper_de$probe)]
}
cat("    Probes with gene symbol:", sum(!is.na(symbol)), "\n")

targets <- union(intersect(TARGET_PROBES, rownames(expr)),
                 names(symbol)[symbol %in% TARGET_GENES])
if (length(targets) == 0) stop("No PON2/SQSTM1 probe found in expression data.")
target_label <- setNames(paste0(symbol[targets], "\n(", targets, ")"), targets)
target_label[is.na(symbol[targets])] <- targets[is.na(symbol[targets])]
cat("    Target probes:", paste(targets, symbol[targets], sep = "=", collapse = ", "), "\n\n")

# ── 5. Sample metadata ───────────────────────────────────────────────────────
cat("[3] Building sample metadata ...\n")
pick_ch1 <- function(pattern) {
  col <- grep(pattern, colnames(pheno), ignore.case = TRUE, value = TRUE)
  col <- col[grepl(":ch1$", col)]
  if (length(col) > 0) as.character(pheno[[col[1]]]) else rep(NA_character_, nrow(pheno))
}
norm_sex <- function(x) case_when(
  grepl("^\\s*(f|female|woman)", x, ignore.case = TRUE) ~ "Female",
  grepl("^\\s*(m|male|man)",     x, ignore.case = TRUE) ~ "Male",
  TRUE ~ NA_character_)

id_text <- apply(pheno[, intersect(c("title", "description", "source_name_ch1"),
                                   colnames(pheno)), drop = FALSE],
                 1, paste, collapse = " ")
id_m <- str_match(id_text, "(?i)(?<![A-Za-z])(CT|PD)[ _.-]?0*([0-9]{1,2})(?![0-9])")
geo_grp_txt <- pick_ch1("disease|diagnosis|group|condition|status|phenotype")

meta <- tibble(
  geo_accession = rownames(pheno),
  title         = as.character(pheno$title),
  sample_id     = ifelse(is.na(id_m[, 1]), NA, paste0(toupper(id_m[, 2]), id_m[, 3])),
  geo_age       = suppressWarnings(as.numeric(str_extract(pick_ch1("^age"), "[0-9.]+"))),
  geo_sex       = norm_sex(pick_ch1("^(sex|gender)")),
  geo_group     = case_when(
    grepl("parkinson|\\bPD\\b", geo_grp_txt, ignore.case = TRUE)        ~ "PD",
    grepl("control|healthy|\\bCT\\b|\\bHC\\b", geo_grp_txt, ignore.case = TRUE) ~ "CT",
    TRUE ~ NA_character_)
) %>%
  left_join(paper_meta, by = "sample_id") %>%
  mutate(
    group = coalesce(geo_group, str_extract(sample_id, "^(CT|PD)")),
    age   = coalesce(geo_age, age),
    sex   = coalesce(geo_sex, sex)
  )

if (!is.null(SAMPLE_META_FILE)) {
  ov <- read_csv(SAMPLE_META_FILE, show_col_types = FALSE)
  meta <- meta %>%
    rows_update(ov %>% select(any_of(c("geo_accession", "group", "sex", "age", "site", "rin"))),
                by = "geo_accession", unmatched = "ignore")
}

meta <- meta %>%
  filter(geo_accession %in% colnames(expr)) %>%
  mutate(group = factor(group, levels = c("CT", "PD")),
         sex   = factor(sex,   levels = c("Female", "Male")))

if (any(is.na(meta$group) | is.na(meta$sex) | is.na(meta$age))) {
  print(meta %>% filter(is.na(group) | is.na(sex) | is.na(age)) %>%
          select(geo_accession, title, sample_id, group, sex, age))
  stop("Group/sex/age missing for the samples above. ",
       "Provide SAMPLE_META_FILE (geo_accession, group, sex, age).")
}
expr <- expr[, meta$geo_accession]

print(meta %>% group_by(group) %>%
        summarise(n = n(), female = sum(sex == "Female"),
                  age_mean = mean(age), age_sd = sd(age),
                  age_min = min(age), age_max = max(age)))
cat(sprintf("    Age difference (Welch t): p = %.2g\n\n",
            t.test(age ~ group, data = meta)$p.value))
write_csv(meta, file.path(OUTPUT_DIR, "sample_metadata.csv"))

# ── 6. Probe filter (as in paper) ────────────────────────────────────────────
expr_all <- expr
if (FILTER_TO_PAPER_PROBES && !is.null(paper_de)) {
  keep <- rownames(expr) %in% paper_de$probe
  cat("[4] Filter: probes in paper table ->", sum(keep), "\n")
} else {
  # Paper: ~18,091 of ~21,448 genes kept (low-expression filter)
  med  <- apply(expr, 1, median)
  keep <- med >= quantile(med, 0.16)
  cat("[4] Filter: row median >= 16th percentile ->", sum(keep), "\n")
}
keep[rownames(expr) %in% targets] <- TRUE
expr <- expr[keep, , drop = FALSE]
cat("\n")

# ── 7. Cell-composition covariates ───────────────────────────────────────────
cat("[5] Cell-composition covariates ...\n")
if (!is.null(CIBERSORT_FILE)) {
  cs <- if (grepl("\\.xlsx?$", CIBERSORT_FILE)) read_excel(CIBERSORT_FILE)
        else read_delim(CIBERSORT_FILE, show_col_types = FALSE)
  names(cs)[1] <- "mixture"
  cs <- cs %>% select(-any_of(c("P-value", "Correlation", "RMSE", "Absolute score (sig.score)")))
  names(cs) <- make.names(names(cs))
  cs$mixture <- as.character(cs$mixture)
  key <- case_when(
    cs$mixture %in% meta$geo_accession ~ cs$mixture,
    cs$mixture %in% meta$sample_id     ~ meta$geo_accession[match(cs$mixture, meta$sample_id)],
    TRUE ~ str_extract(cs$mixture, "GSM[0-9]+"))
  cell_scores <- cs %>% mutate(geo_accession = key) %>%
    select(-mixture) %>% filter(!is.na(geo_accession))
  COMP_COVARS <- c("Neutrophils", "T.cells.CD4.naive")
  cat("    Using CIBERSORTx proportions from", CIBERSORT_FILE, "\n")
} else {
  # Marker-gene scores: mean z-score of marker probes (relative abundance
  # proxy, not absolute proportion). T_naive is a proxy for naive CD4 T.
  MARKERS <- list(
    Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "CXCR1", "FPR1", "CMTM2", "MME", "ALPL"),
    T_naive    = c("CCR7", "LEF1", "SELL", "TCF7", "MAL", "NOSIP"),
    T_CD4      = c("CD4", "IL7R", "CD40LG"),
    T_CD8      = c("CD8A", "CD8B"),
    Monocyte   = c("CD14", "LYZ", "FCN1", "CSF1R", "VCAN"),
    NK         = c("GNLY", "NKG7", "KLRD1", "KLRF1", "PRF1"),
    B_cell     = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
    Eosinophil = c("CLC", "IL5RA", "CCR3", "SIGLEC8"),
    Platelet   = c("PPBP", "PF4", "ITGA2B", "GP9")
  )
  z <- t(scale(t(expr_all)))
  cell_scores <- map_dfc(MARKERS, function(g) {
    p <- names(symbol)[symbol %in% g]
    if (length(p) == 0) return(rep(NA_real_, ncol(z)))
    colMeans(z[p, , drop = FALSE])
  }) %>% mutate(geo_accession = colnames(z))
  map_int(MARKERS, ~ sum(symbol %in% .x)) %>%
    { cat("    Marker probes found:", paste(names(.), ., sep = "=", collapse = ", "), "\n") }
  cell_scores <- cell_scores %>% select(where(~ !all(is.na(.x))))
  COMP_COVARS <- c("Neutrophil", "T_naive")
  cat("    Using marker-gene scores (set CIBERSORT_FILE to use CIBERSORTx)\n")
}
COMP_COVARS <- intersect(COMP_COVARS, names(cell_scores))
cell_types  <- setdiff(names(cell_scores), "geo_accession")

meta <- meta %>% left_join(cell_scores, by = "geo_accession")

# PCs of the whole composition (summarises all cell types in 2 covariates)
comp_mat <- as.matrix(meta[, cell_types])
comp_mat <- comp_mat[, apply(comp_mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
comp_pca <- prcomp(comp_mat, scale. = TRUE)
meta$compPC1 <- comp_pca$x[, 1]
meta$compPC2 <- comp_pca$x[, 2]
cat(sprintf("    Composition PC1/PC2 explain %.1f%% / %.1f%%\n",
            100 * summary(comp_pca)$importance[2, 1],
            100 * summary(comp_pca)$importance[2, 2]))

comp_test <- map_dfr(cell_types, function(ct) {
  tt <- t.test(meta[[ct]] ~ meta$group)
  tibble(cell_type = ct,
         mean_CT = tt$estimate[1], mean_PD = tt$estimate[2],
         p_value = tt$p.value,
         r_age   = cor(meta[[ct]], meta$age))
})
print(comp_test, n = Inf)
write_csv(comp_test, file.path(OUTPUT_DIR, "cell_composition_PD_vs_CT.csv"))
cat("\n")

# ── 8. limma models ───────────────────────────────────────────────────────────
cov_str <- paste(COMP_COVARS, collapse = " + ")
MODELS <- list(
  M0_paper_sex        = "~ group + sex",
  M1_sex_age          = "~ group + sex + age",
  M2_sex_cells        = paste("~ group + sex +", cov_str),
  M3_sex_age_cells    = paste("~ group + sex + age +", cov_str),
  M4_sex_age_compPC   = "~ group + sex + age + compPC1 + compPC2"
)
if (all(!is.na(meta$site)) && length(unique(meta$site)) > 1 && all(!is.na(meta$rin))) {
  MODELS$M5_M3_site_RIN <- paste(MODELS$M3_sex_age_cells, "+ site + rin")
}

run_limma <- function(form, m, e, label) {
  design <- model.matrix(as.formula(form), data = m)
  colnames(design) <- make.names(colnames(design))
  fit <- eBayes(lmFit(e[, m$geo_accession], design))

  # Collinearity of group with the other covariates
  others <- setdiff(colnames(design), c("X.Intercept.", "groupPD"))
  vif <- if (length(others) > 0) {
    r2 <- summary(lm(design[, "groupPD"] ~ design[, others]))$r.squared
    1 / (1 - r2)
  } else 1

  tt <- topTable(fit, coef = "groupPD", number = Inf, sort.by = "none", confint = TRUE)
  res <- tt[targets, ] %>%
    rownames_to_column("probe") %>%
    transmute(analysis = label, formula = form, probe, gene = symbol[probe],
              n_CT = sum(m$group == "CT"), n_PD = sum(m$group == "PD"),
              logFC, CI.L, CI.R, t, P.Value, adj.P.Val, VIF_group = vif) %>%
    as_tibble()

  cov_eff <- map_dfr(colnames(design)[-1], function(cf) tibble(
    analysis = label, probe = targets, gene = symbol[targets], term = cf,
    estimate = fit$coefficients[targets, cf],
    p_value  = fit$p.value[targets, cf]))

  list(res = res, cov = cov_eff, fit = fit, design = design,
       genome = tt %>% rownames_to_column("probe") %>% select(probe, logFC, P.Value))
}

cat("[6] Running models (full cohort) ...\n")
fits <- imap(MODELS, ~ run_limma(.x, meta, expr, .y))

# ── 8a. Reproduction check vs paper ──────────────────────────────────────────
if (!is.null(paper_de)) {
  chk <- fits$M0_paper_sex$genome %>%
    inner_join(paper_de, by = "probe", suffix = c("_ours", "_paper"))
  cat(sprintf("    M0 vs paper genome-wide logFC: r = %.3f (n = %d probes)\n",
              cor(chk$logFC_ours, chk$logFC_paper), nrow(chk)))
  cat("    (r close to 1 = labels and normalization match the paper)\n")
  print(chk %>% filter(probe %in% targets) %>%
          select(probe, GENE_NAME, logFC_ours, logFC_paper, P.Value_ours, P.Value_paper))
}

# ── 9. Sensitivity analyses ──────────────────────────────────────────────────
cat("\n[7] Sensitivity analyses ...\n")
sens <- list()

# 9a. Age-overlap subset (both groups within the shared age range)
age_lo <- max(tapply(meta$age, meta$group, min))
age_hi <- min(tapply(meta$age, meta$group, max))
m_ov   <- meta %>% filter(age >= age_lo, age <= age_hi)
cat(sprintf("    Age-overlap subset [%g, %g]: CT = %d, PD = %d\n",
            age_lo, age_hi, sum(m_ov$group == "CT"), sum(m_ov$group == "PD")))
if (min(table(m_ov$group)) >= 3) {
  sens$S1_ageOverlap_sex     <- run_limma("~ group + sex",       m_ov, expr, "S1_ageOverlap_sex")
  sens$S1_ageOverlap_sex_age <- run_limma("~ group + sex + age", m_ov, expr, "S1_ageOverlap_sex_age")
}

# 9b. Remove extreme age (PD7, 106 y) — high-leverage point for the age term
m_noout <- meta %>% filter(age < 100)
if (nrow(m_noout) < nrow(meta)) {
  sens$S2_noAgeOutlier_M3 <- run_limma(MODELS$M3_sex_age_cells, m_noout, expr, "S2_noAgeOutlier_M3")
}

# 9c. Age slope within each group (age effect without disease contrast)
age_within <- map_dfr(c("CT", "PD"), function(g) {
  mg <- meta %>% filter(group == g)
  d  <- model.matrix(~ sex + age, data = mg)
  f  <- eBayes(lmFit(expr[, mg$geo_accession], d))
  tibble(group = g, n = nrow(mg), probe = targets, gene = symbol[targets],
         age_slope_per_10y = 10 * f$coefficients[targets, "age"],
         p_value = f$p.value[targets, "age"])
})
cat("    Age slope (log2 per 10 years) within groups:\n")
print(age_within)
write_csv(age_within, file.path(OUTPUT_DIR, "age_slope_within_group.csv"))

# ── 10. Collect results ──────────────────────────────────────────────────────
all_fits   <- c(fits, sens)
results    <- map_dfr(all_fits, "res")
cov_effect <- map_dfr(all_fits, "cov")
write_csv(results,    file.path(OUTPUT_DIR, "PON2_SQSTM1_PD_effect_by_model.csv"))
write_csv(cov_effect, file.path(OUTPUT_DIR, "PON2_SQSTM1_covariate_effects.csv"))

cat("\n[8] PD vs CT effect across models:\n")
print(results %>%
        transmute(analysis, gene, probe, n_CT, n_PD,
                  logFC = round(logFC, 3),
                  CI = sprintf("[%.3f, %.3f]", CI.L, CI.R),
                  P.Value = signif(P.Value, 3), adj.P.Val = signif(adj.P.Val, 3),
                  VIF_group = round(VIF_group, 2)),
      n = Inf, width = Inf)

attn <- results %>%
  filter(analysis %in% c("M0_paper_sex", "M1_sex_age", "M2_sex_cells", "M3_sex_age_cells")) %>%
  select(analysis, probe, gene, logFC) %>%
  pivot_wider(names_from = analysis, values_from = logFC) %>%
  mutate(pct_change_age        = 100 * (M1_sex_age       / M0_paper_sex - 1),
         pct_change_cells      = 100 * (M2_sex_cells     / M0_paper_sex - 1),
         pct_change_age_cells  = 100 * (M3_sex_age_cells / M0_paper_sex - 1))
cat("\n    % change of logFC relative to paper model (negative = attenuated):\n")
print(attn %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)
write_csv(attn, file.path(OUTPUT_DIR, "logFC_attenuation.csv"))

# ── 11. Plots ─────────────────────────────────────────────────────────────────
cat("\n[9] Plotting ...\n")
grp_cols <- c(CT = "#4C78A8", PD = "#E45756")

long_expr <- as_tibble(t(expr[targets, , drop = FALSE]), rownames = "geo_accession") %>%
  pivot_longer(-geo_accession, names_to = "probe", values_to = "expr") %>%
  left_join(meta, by = "geo_accession") %>%
  mutate(label = target_label[probe])

# 11a. Forest plot of PD effect across models
paper_ref <- if (!is.null(paper_de)) {
  paper_de %>% filter(probe %in% targets) %>% transmute(probe, paper_logFC = logFC)
} else tibble(probe = character(), paper_logFC = numeric())

p_forest <- results %>%
  mutate(label = target_label[probe],
         analysis = factor(analysis, levels = rev(unique(analysis))),
         sig = ifelse(P.Value < 0.05, "P < 0.05", "n.s.")) %>%
  ggplot(aes(x = logFC, y = analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(data = paper_ref %>% mutate(label = target_label[probe]),
             aes(xintercept = paper_logFC), colour = "#E45756", alpha = 0.5) +
  geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), height = 0.25) +
  geom_point(aes(fill = sig), shape = 21, size = 2.8) +
  scale_fill_manual(values = c("P < 0.05" = "black", "n.s." = "white"), name = NULL) +
  facet_wrap(~ label, scales = "free_x") +
  labs(title = "PD vs CT effect (log2FC, 95% CI) by adjustment model",
       subtitle = "Red line = paper estimate (~ group + sex)",
       x = "log2 fold change (PD - CT)", y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "01_forest_PD_effect_by_model.png"),
       p_forest, width = 4 + 3.5 * length(targets), height = 5.5, dpi = 300)

# 11b. Expression vs age
p_age <- ggplot(long_expr, aes(age, expr, colour = group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, formula = y ~ x,
              colour = "grey30", linetype = "dashed") +
  scale_colour_manual(values = grp_cols) +
  facet_wrap(~ label, scales = "free_y") +
  labs(title = "Expression vs age", subtitle = "Dashed = all samples",
       x = "Age (years)", y = "log2 expression (RMA)") +
  theme_bw(base_size = 11)
ggsave(file.path(OUTPUT_DIR, "02_expression_vs_age.png"),
       p_age, width = 3 + 3.5 * length(targets), height = 4.5, dpi = 300)

# 11c. Raw vs covariate-adjusted expression (M3: sex + age + cells removed)
adjusted_expr <- function(f, m) {
  d    <- f$design
  covs <- setdiff(colnames(d), c("X.Intercept.", "groupPD"))
  b    <- f$fit$coefficients[targets, covs, drop = FALSE]
  xc   <- sweep(d[, covs, drop = FALSE], 2, colMeans(d[, covs, drop = FALSE]))
  expr[targets, m$geo_accession, drop = FALSE] - b %*% t(xc)
}
adj_m3 <- adjusted_expr(fits$M3_sex_age_cells, meta)
long_adj <- as_tibble(t(adj_m3), rownames = "geo_accession") %>%
  pivot_longer(-geo_accession, names_to = "probe", values_to = "expr") %>%
  left_join(meta %>% select(geo_accession, group), by = "geo_accession") %>%
  mutate(type = "Adjusted (M3: sex + age + cells)")

box_df <- bind_rows(long_expr %>% select(geo_accession, probe, expr, group) %>%
                      mutate(type = "Raw"),
                    long_adj) %>%
  mutate(label = target_label[probe],
         type  = factor(type, levels = c("Raw", "Adjusted (M3: sex + age + cells)")))
p_lab <- results %>%
  filter(analysis %in% c("M0_paper_sex", "M3_sex_age_cells")) %>%
  mutate(type = ifelse(analysis == "M0_paper_sex", "Raw", "Adjusted (M3: sex + age + cells)"),
         type = factor(type, levels = levels(box_df$type)),
         label = target_label[probe],
         txt = sprintf("logFC = %.2f\nP = %.2g", logFC, P.Value))

p_box <- ggplot(box_df, aes(group, expr, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.8) +
  geom_text(data = p_lab, aes(x = 1.5, y = Inf, label = txt), inherit.aes = FALSE,
            vjust = 1.2, size = 3) +
  scale_fill_manual(values = grp_cols) +
  facet_grid(type ~ label, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(title = "PON2 / SQSTM1: raw (paper model P) vs adjusted",
       x = NULL, y = "log2 expression") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(OUTPUT_DIR, "03_boxplot_raw_vs_adjusted.png"),
       p_box, width = 2 + 3 * length(targets), height = 7, dpi = 300)

# 11d. Cell composition by group + correlation with targets
p_comp <- meta %>%
  select(geo_accession, group, all_of(cell_types)) %>%
  pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "score") %>%
  ggplot(aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1) +
  scale_fill_manual(values = grp_cols) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
  labs(title = "Cell-composition covariates by group",
       y = if (is.null(CIBERSORT_FILE)) "Marker score (mean z)" else "CIBERSORTx fraction",
       x = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

cor_df <- expand_grid(probe = targets, covariate = c("age", cell_types)) %>%
  mutate(r = map2_dbl(probe, covariate, ~ cor(expr[.x, meta$geo_accession], meta[[.y]])),
         label = gsub("\n", " ", target_label[probe]))
p_cor <- ggplot(cor_df, aes(covariate, label, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#4C78A8", high = "#E45756", limits = c(-1, 1)) +
  labs(title = "Pearson r: target expression vs age / cell covariates", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "04_cell_composition_and_correlation.png"),
       p_comp / p_cor + plot_layout(heights = c(2, 1)),
       width = 11, height = 9, dpi = 300)
write_csv(cor_df %>% select(-label), file.path(OUTPUT_DIR, "correlation_target_vs_covariates.csv"))

# ── 12. Interpretation guide ─────────────────────────────────────────────────
cat("\n=== Interpretation guide ===\n")
cat("* M0 should reproduce the paper (check r and target logFC above).\n")
cat("* If logFC shrinks toward 0 in M1 and the age slope is similar in CT-only\n",
    "  and PD-only fits -> the paper's signal is largely an age effect.\n")
cat("* If logFC shrinks in M2 and the target correlates with Neutrophil /\n",
    "  T_naive scores -> the signal likely reflects cell composition.\n")
cat("* VIF_group > ~2 means group is collinear with the covariates: CIs widen,\n",
    "  so a non-significant M3 is 'not resolvable' rather than 'no effect'.\n",
    "  Compare estimates (forest plot), not only P-values.\n")
cat("* S1 (age-overlap subset) is the design-based check: no extrapolation\n",
    "  of the age slope, but small n.\n")
cat("\nOutputs written to:", normalizePath(OUTPUT_DIR), "\n")
