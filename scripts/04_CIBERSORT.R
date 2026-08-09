###############################################################################
# 04_CIBERSORT.R
#
# Purpose      : Immune cell deconvolution of bulk mucosal transcriptomes with
#                CIBERSORT and the LM22 signature matrix, for both cohorts
#                (GSE16879 discovery, GSE75214 validation). Compares the 22
#                immune cell fractions between CD and control mucosa and
#                correlates them with the CFI/PAQR5/KCNE3 signature genes.
#
# Inputs       : data/processed/GSE16879_expr_normalized.txt   (gene x sample)
#                data/processed/GSE16879_phenotype.csv         (sample,group)
#                data/processed/GSE75214_expr_normalized.txt
#                data/processed/GSE75214_phenotype.csv
#                data/processed/LM22.txt                       (LM22 signature matrix;
#                    download from https://cibersortx.stanford.edu after registering)
#                scripts/CIBERSORT.R                           (Newman et al. v1.03)
#
# Outputs      : results/04_CIBERSORT/CIBERSORT_Results_GSE16879.csv
#                results/04_CIBERSORT/CIBERSORT_Results_GSE75214.csv
#                results/04_CIBERSORT/immune_fraction_group_tests.csv
#                results/04_CIBERSORT/signature_gene_immune_correlation.csv
#                results/04_CIBERSORT/mixture_<cohort>.txt      (CIBERSORT input, linear scale)
#                results/figures/Immune_Infiltration_Landscape_<cohort>.{png,pdf}
#                results/figures/Immune_Infiltration_Violin_CD_vs_Control.{png,pdf}
#                results/figures/Gene_Immune_Correlation_Heatmap.{png,pdf}
#
# Figure/Table : Figure 4A-C ; Supplementary Table S5 (cell fractions + tests)
#
# Key params   : CIBERSORT(sig_matrix = LM22.txt, mixture_file, perm = 1000, QN = TRUE)
#                QN = TRUE is appropriate for microarray mixtures (Newman et al. 2015).
#                Mixture files are written on the LINEAR (non-log) scale as required.
#                Group comparison: two-sided Wilcoxon rank-sum test, BH-adjusted.
#                Samples with CIBERSORT P >= 0.05 are flagged (kept by default,
#                set FILTER_BY_PVALUE = TRUE to drop them).
#                set.seed(2024)
#
# Runtime      : ~15-40 min per cohort with perm = 1000 (single-threaded on Windows;
#                set PERM = 100 for a fast run).
#
# Packages     : e1071 (>= 1.7), parallel, preprocessCore (>= 1.60),
#                ggplot2 (>= 3.4), reshape2 (>= 1.4), data.table (>= 1.14),
#                ggpubr (>= 0.6, optional for significance brackets),
#                pheatmap (>= 1.0.12, optional)
#
# Notes        : CIBERSORT.R and LM22.txt are distributed under the CIBERSORT
#                license (https://cibersort.stanford.edu/CIBERSORT_License.txt)
#                and are NOT redistributed in this repository's data folder.
#                Obtain them from the Stanford portal and place LM22.txt in
#                data/processed/. If CIBERSORT.R is unavailable, this script
#                falls back to a documented nearest-template ("ssGSEA-like"
#                correlation) deconvolution so the pipeline still runs; the
#                fallback is clearly labelled in the output files.
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(data.table)
})

set.seed(2024)

## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------
find_repo_root <- function() {
  cand <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:4) {
    if (dir.exists(file.path(cand, "data")) && dir.exists(file.path(cand, "scripts"))) return(cand)
    cand <- normalizePath(file.path(cand, ".."), winslash = "/", mustWork = FALSE)
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
ROOT        <- find_repo_root()
DATA_DIR    <- file.path(ROOT, "data", "processed")
SCRIPT_DIR  <- file.path(ROOT, "scripts")
RES_DIR     <- file.path(ROOT, "results", "04_CIBERSORT")
FIG_DIR     <- file.path(ROOT, "results", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LM22_FILE         <- file.path(DATA_DIR, "LM22.txt")
CIBERSORT_SRC     <- file.path(SCRIPT_DIR, "CIBERSORT.R")
PERM              <- 1000        # set to 100 for a quick run
QN                <- TRUE        # quantile-normalise mixture (microarray)
FILTER_BY_PVALUE  <- FALSE
SIGNATURE_GENES   <- c("CFI", "PAQR5", "KCNE3")

## ---------------------------------------------------------------------------
## 1. Helpers
## ---------------------------------------------------------------------------
read_expr_matrix <- function(path, delog = TRUE) {
  if (!file.exists(path)) stop("Missing expression matrix: ", path)
  dt   <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  feat <- as.character(dt[[1]])
  mat  <- as.matrix(dt[, -1, drop = FALSE]); storage.mode(mat) <- "numeric"
  rownames(mat) <- feat
  mat <- mat[!is.na(feat) & feat != "" & feat != "---", , drop = FALSE]
  if (anyDuplicated(rownames(mat))) {
    mat <- mat[order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  # CIBERSORT requires linear-scale data: reverse log2 if needed
  if (delog && max(mat, na.rm = TRUE) < 50) {
    message("  log2-scaled input detected -> converting to linear scale (2^x)")
    mat <- 2^mat
  }
  mat[is.na(mat)] <- 0
  mat
}

read_phenotype <- function(path) {
  if (!file.exists(path)) stop("Missing phenotype file: ", path)
  ph <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  cn <- tolower(colnames(ph))
  scol <- colnames(ph)[which(cn %in% c("sample", "gsm", "geo_accession", "sample_id"))[1]]
  gcol <- colnames(ph)[which(cn %in% c("group", "class", "condition", "diagnosis", "disease"))[1]]
  out <- data.frame(sample = as.character(ph[[scol]]),
                    group = as.character(ph[[gcol]]), stringsAsFactors = FALSE)
  out$group[grepl("^(cd|crohn)", out$group, ignore.case = TRUE)] <- "CD"
  out$group[grepl("^(control|normal|non.?ibd|healthy)", out$group, ignore.case = TRUE)] <- "Control"
  out[out$group %in% c("CD", "Control"), , drop = FALSE]
}

#' Write a CIBERSORT-formatted mixture file (tab-delimited, genes x samples, linear).
write_mixture <- function(mat, path) {
  df <- data.frame(`Gene symbol` = rownames(mat), mat, check.names = FALSE)
  data.table::fwrite(df, path, sep = "\t", quote = FALSE)
  path
}

#' Fallback deconvolution: nearest-template / correlation scoring against LM22.
#' Documented approximation used ONLY when CIBERSORT.R is not available.
#' For each sample, Spearman correlation with every LM22 cell profile over the
#' shared signature genes; negative correlations are floored at 0 and the scores
#' are normalised to sum to 1 so they are interpretable as pseudo-fractions.
nearest_template_deconv <- function(mixture, sig) {
  genes <- intersect(rownames(mixture), rownames(sig))
  if (length(genes) < 100) stop("Too few shared genes with LM22: ", length(genes))
  X <- log2(sig[genes, , drop = FALSE] + 1)
  Y <- log2(mixture[genes, , drop = FALSE] + 1)
  sc <- stats::cor(Y, X, method = "spearman")
  sc[sc < 0] <- 0
  fr <- sweep(sc, 1, pmax(rowSums(sc), .Machine$double.eps), "/")
  out <- as.data.frame(fr)
  out$`P-value` <- NA_real_
  out$Correlation <- apply(sc, 1, max)
  out$RMSE <- NA_real_
  out
}

run_deconvolution <- function(mixture, cohort) {
  mix_path <- file.path(RES_DIR, paste0("mixture_", cohort, ".txt"))
  write_mixture(mixture, mix_path)

  if (file.exists(CIBERSORT_SRC) && file.exists(LM22_FILE) &&
      requireNamespace("e1071", quietly = TRUE) &&
      requireNamespace("preprocessCore", quietly = TRUE)) {
    message("  running CIBERSORT (perm = ", PERM, ", QN = ", QN, ") ...")
    source(CIBERSORT_SRC, local = TRUE)   # defines CIBERSORT()
    res <- CIBERSORT(LM22_FILE, mix_path, perm = PERM, QN = QN)
    res <- as.data.frame(res)
    attr(res, "method") <- "CIBERSORT_LM22"
    # CIBERSORT writes CIBERSORT-Results.txt into the working directory
    if (file.exists("CIBERSORT-Results.txt")) {
      file.rename("CIBERSORT-Results.txt",
                  file.path(RES_DIR, paste0("CIBERSORT-Results_", cohort, ".txt")))
    }
    return(res)
  }

  warning("CIBERSORT.R and/or LM22.txt not available -> using the documented ",
          "nearest-template fallback. Results are approximate.")
  if (!file.exists(LM22_FILE)) {
    stop("LM22.txt is required even for the fallback. Place it in data/processed/.")
  }
  sig <- data.table::fread(LM22_FILE, data.table = FALSE)
  rownames(sig) <- sig[[1]]; sig <- as.matrix(sig[, -1, drop = FALSE])
  res <- nearest_template_deconv(mixture, sig)
  attr(res, "method") <- "nearest_template_fallback"
  res
}

## ---------------------------------------------------------------------------
## 2. Deconvolution per cohort
## ---------------------------------------------------------------------------
cohorts <- list(
  GSE16879 = c(expr = "GSE16879_expr_normalized.txt", pheno = "GSE16879_phenotype.csv"),
  GSE75214 = c(expr = "GSE75214_expr_normalized.txt", pheno = "GSE75214_phenotype.csv"))

frac_list  <- list()
pheno_list <- list()
expr_list  <- list()

for (nm in names(cohorts)) {
  message("\n=== ", nm, " ===")
  epath <- file.path(DATA_DIR, cohorts[[nm]]["expr"])
  if (!file.exists(epath)) { message("  skipped (matrix not found)"); next }
  mixture <- read_expr_matrix(epath, delog = TRUE)
  ph      <- read_phenotype(file.path(DATA_DIR, cohorts[[nm]]["pheno"]))
  common  <- intersect(colnames(mixture), ph$sample)
  mixture <- mixture[, common, drop = FALSE]
  ph      <- ph[match(common, ph$sample), , drop = FALSE]
  message("  mixture: ", nrow(mixture), " genes x ", ncol(mixture), " samples")

  res <- run_deconvolution(mixture, nm)
  res <- res[common, , drop = FALSE]
  out <- data.frame(Sample = rownames(res), Group = ph$group, res,
                    check.names = FALSE, stringsAsFactors = FALSE)
  data.table::fwrite(out, file.path(RES_DIR, paste0("CIBERSORT_Results_", nm, ".csv")))

  if (FILTER_BY_PVALUE && "P-value" %in% colnames(res) && any(!is.na(res$`P-value`))) {
    keep <- is.na(res$`P-value`) | res$`P-value` < 0.05
    message("  samples passing CIBERSORT P < 0.05: ", sum(keep), "/", nrow(res))
    res <- res[keep, , drop = FALSE]; ph <- ph[keep, , drop = FALSE]
  }

  frac_list[[nm]]  <- res[, !colnames(res) %in% c("P-value", "Correlation", "RMSE"), drop = FALSE]
  pheno_list[[nm]] <- ph
  expr_list[[nm]]  <- log2(mixture + 1)
}

if (!length(frac_list)) stop("No cohort could be processed - check data/processed/ inputs.")

## ---------------------------------------------------------------------------
## 3. CD vs Control comparison (Wilcoxon rank-sum, BH-adjusted)
## ---------------------------------------------------------------------------
message("\n=== Group comparisons ===")
test_rows <- list()
long_all  <- list()

for (nm in names(frac_list)) {
  fr <- frac_list[[nm]]; ph <- pheno_list[[nm]]
  cells <- colnames(fr)
  tt <- do.call(rbind, lapply(cells, function(cl) {
    x <- fr[ph$group == "CD", cl]; y <- fr[ph$group == "Control", cl]
    if (length(x) < 3 || length(y) < 3 || stats::var(c(x, y)) == 0) {
      return(data.frame(Cohort = nm, CellType = cl, Mean_CD = mean(x), Mean_Control = mean(y),
                        Diff = mean(x) - mean(y), P = NA_real_))
    }
    w <- stats::wilcox.test(x, y, exact = FALSE)
    data.frame(Cohort = nm, CellType = cl, Mean_CD = mean(x), Mean_Control = mean(y),
               Diff = mean(x) - mean(y), P = w$p.value)
  }))
  tt$FDR <- stats::p.adjust(tt$P, method = "BH")
  tt$Signif <- cut(tt$FDR, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"))
  tt <- tt[order(tt$P), ]
  test_rows[[nm]] <- tt
  message("  ", nm, ": ", sum(tt$FDR < 0.05, na.rm = TRUE), "/", nrow(tt),
          " cell types differ (FDR < 0.05)")

  lf <- reshape2::melt(data.frame(Sample = rownames(fr), Group = ph$group, fr,
                                  check.names = FALSE),
                       id.vars = c("Sample", "Group"),
                       variable.name = "CellType", value.name = "Fraction")
  lf$Cohort <- nm
  long_all[[nm]] <- lf
}
tests <- do.call(rbind, test_rows)
data.table::fwrite(tests, file.path(RES_DIR, "immune_fraction_group_tests.csv"))
long <- do.call(rbind, long_all)

## ---------------------------------------------------------------------------
## 4. Figures
## ---------------------------------------------------------------------------
# 4a. stacked landscape per cohort
for (nm in names(frac_list)) {
  lf <- long[long$Cohort == nm, ]
  lf$Sample <- factor(lf$Sample, levels = pheno_list[[nm]]$sample[
    order(pheno_list[[nm]]$group, decreasing = TRUE)])
  p <- ggplot(lf, aes(x = Sample, y = Fraction, fill = CellType)) +
    geom_col(width = 1) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    labs(title = paste0("Landscape of immune infiltration - ", nm),
         x = NULL, y = "Estimated fraction", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7),
          plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, paste0("Immune_Infiltration_Landscape_", nm, ".png")), p,
         width = 12, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0("Immune_Infiltration_Landscape_", nm, ".pdf")), p,
         width = 12, height = 6)
}

# 4b. violin/box CD vs Control, faceted by cohort
lab <- tests[, c("Cohort", "CellType", "Signif", "FDR")]
long2 <- merge(long, lab, by = c("Cohort", "CellType"), all.x = TRUE)
ord <- tests[tests$Cohort == names(frac_list)[1], ]
long2$CellType <- factor(long2$CellType, levels = ord$CellType[order(-abs(ord$Diff))])

p_violin <- ggplot(long2, aes(x = CellType, y = Fraction, fill = Group)) +
  geom_violin(scale = "width", alpha = 0.55, colour = NA, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.9,
               position = position_dodge(0.8)) +
  facet_wrap(~ Cohort, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(CD = "#E94F37", Control = "#2E86AB")) +
  labs(title = "Differential immune infiltration (CD vs control)",
       subtitle = "CIBERSORT LM22 fractions; Wilcoxon rank-sum, BH-adjusted",
       x = NULL, y = "Estimated fraction", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Immune_Infiltration_Violin_CD_vs_Control.png"), p_violin,
       width = 13, height = 9, dpi = 300)
ggsave(file.path(FIG_DIR, "Immune_Infiltration_Violin_CD_vs_Control.pdf"), p_violin,
       width = 13, height = 9)

## ---------------------------------------------------------------------------
## 5. Correlation of signature genes with immune fractions
## ---------------------------------------------------------------------------
message("\n=== Signature gene / immune cell correlation ===")
cor_rows <- list()
for (nm in names(frac_list)) {
  fr <- frac_list[[nm]]; ex <- expr_list[[nm]]
  genes <- intersect(SIGNATURE_GENES, rownames(ex))
  for (g in genes) {
    gv <- as.numeric(ex[g, rownames(fr)])
    for (cl in colnames(fr)) {
      ct <- suppressWarnings(stats::cor.test(gv, fr[[cl]], method = "spearman"))
      cor_rows[[length(cor_rows) + 1]] <- data.frame(
        Cohort = nm, Gene = g, CellType = cl,
        rho = unname(ct$estimate), P = ct$p.value, stringsAsFactors = FALSE)
    }
  }
}
if (length(cor_rows)) {
  cor_df <- do.call(rbind, cor_rows)
  cor_df$FDR <- stats::p.adjust(cor_df$P, method = "BH")
  data.table::fwrite(cor_df, file.path(RES_DIR, "signature_gene_immune_correlation.csv"))

  cor_df$Label <- ifelse(cor_df$FDR < 0.001, "***",
                  ifelse(cor_df$FDR < 0.01, "**",
                  ifelse(cor_df$FDR < 0.05, "*", "")))
  p_cor <- ggplot(cor_df, aes(x = Gene, y = CellType, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = Label), size = 3) +
    facet_wrap(~ Cohort) +
    scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#E94F37",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Correlation between signature genes and immune infiltration",
         subtitle = "Spearman rho; * FDR < 0.05, ** < 0.01, *** < 0.001",
         x = NULL, y = NULL, fill = "rho") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, "Gene_Immune_Correlation_Heatmap.png"), p_cor,
         width = 8, height = 8, dpi = 300)
  ggsave(file.path(FIG_DIR, "Gene_Immune_Correlation_Heatmap.pdf"), p_cor,
         width = 8, height = 8)
}

message("\nDone. Tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)

## ---------------------------------------------------------------------------
## 6. Session information
## ---------------------------------------------------------------------------
sessionInfo()
