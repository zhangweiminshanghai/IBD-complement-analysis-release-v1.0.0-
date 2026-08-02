###############################################################################
# 08_bootstrap_AUC.R
#
# Purpose      : Internal validation of the three-gene diagnostic signature
#                (CFI, PAQR5, KCNE3) by non-parametric bootstrap resampling.
#                For each of 2,000 bootstrap replicates the logistic model is
#                refitted on the resampled data and evaluated on the out-of-bag
#                (OOB) samples, giving an optimism-corrected estimate of
#                discrimination:
#                    mean bootstrap AUC = 0.711 (95% CI 0.565 - 0.844)
#                The apparent (in-sample) AUC and the .632+ style optimism
#                correction are reported alongside for transparency.
#
# Inputs       : data/processed/GSE16879_expr_normalized.txt  (discovery, linear -> log2)
#                data/processed/GSE16879_phenotype.csv        (sample,group)
#                data/processed/GSE75214_expr_normalized.txt  (optional, for the
#                    external bootstrap distribution)
#                data/processed/GSE75214_phenotype.csv        (optional)
#                results/03_ML/signature_genes.csv            (optional; defaults to
#                    CFI/PAQR5/KCNE3 when absent)
#
# Outputs      : results/08_bootstrap/bootstrap_auc_replicates.csv
#                results/08_bootstrap/bootstrap_auc_summary.csv
#                results/08_bootstrap/bootstrap_ci_table.csv
#                results/figures/Bootstrap_AUC_distribution.{png,pdf}
#                results/figures/Bootstrap_AUC_forest.{png,pdf}
#
# Figure/Table : Figure 3E (bootstrap distribution) ; Table 2 (validated AUC)
#
# Key params   : N_BOOT      = 2000 resamples (with replacement, stratified by group)
#                MODEL       = logistic regression, group ~ CFI + PAQR5 + KCNE3
#                EVALUATION  = out-of-bag samples of each replicate
#                CI          = percentile bootstrap 95% CI (2.5th, 97.5th)
#                set.seed(2024)
#
# Runtime      : ~1-3 min (2,000 replicates on <= 200 samples)
#
# Packages     : pROC (>= 1.18), ggplot2 (>= 3.4), data.table (>= 1.14)
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
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
ROOT     <- find_repo_root()
DATA_DIR <- file.path(ROOT, "data", "processed")
RES_DIR  <- file.path(ROOT, "results", "08_bootstrap")
FIG_DIR  <- file.path(ROOT, "results", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

N_BOOT          <- 2000
CI_LEVEL        <- 0.95
DEFAULT_GENES   <- c("CFI", "PAQR5", "KCNE3")
REF_MEAN_AUC    <- 0.711     # published bootstrap-validated AUC
REF_CI          <- c(0.565, 0.844)

## ---------------------------------------------------------------------------
## 1. Data loading (same conventions as 01/03)
## ---------------------------------------------------------------------------
read_expr_matrix <- function(path) {
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
  if (max(mat, na.rm = TRUE) > 50) { mat[mat < 0] <- 0; mat <- log2(mat + 1) }
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

build_df <- function(expr_file, pheno_file, genes) {
  e  <- read_expr_matrix(file.path(DATA_DIR, expr_file))
  ph <- read_phenotype(file.path(DATA_DIR, pheno_file))
  common <- intersect(colnames(e), ph$sample)
  e  <- e[, common, drop = FALSE]; ph <- ph[match(common, ph$sample), , drop = FALSE]
  miss <- setdiff(genes, rownames(e))
  if (length(miss)) stop("Signature gene(s) missing: ", paste(miss, collapse = ", "))
  d <- as.data.frame(t(e[genes, , drop = FALSE])); colnames(d) <- genes
  d <- as.data.frame(scale(d))
  d$group <- factor(ph$group, levels = c("Control", "CD"))
  d$sample <- ph$sample
  d
}

sig_file <- file.path(ROOT, "results", "03_ML", "signature_genes.csv")
GENES <- if (file.exists(sig_file)) {
  as.character(data.table::fread(sig_file, data.table = FALSE)$Gene)
} else DEFAULT_GENES
message("Signature: ", paste(GENES, collapse = " + "))

dd <- build_df("GSE16879_expr_normalized.txt", "GSE16879_phenotype.csv", GENES)
message("Discovery cohort: ", nrow(dd), " samples (CD = ", sum(dd$group == "CD"),
        ", Control = ", sum(dd$group == "Control"), ")")

## ---------------------------------------------------------------------------
## 2. Apparent (in-sample) performance
## ---------------------------------------------------------------------------
fml  <- stats::as.formula(paste("group ~", paste(GENES, collapse = " + ")))
full <- stats::glm(fml, data = dd, family = stats::binomial())
roc_app <- pROC::roc(dd$group, stats::fitted(full),
                     levels = c("Control", "CD"), direction = "<", quiet = TRUE)
auc_app <- as.numeric(pROC::auc(roc_app))
message("Apparent AUC (no correction): ", round(auc_app, 3))

## ---------------------------------------------------------------------------
## 3. Bootstrap: 2,000 stratified resamples, evaluated out-of-bag
## ---------------------------------------------------------------------------
safe_auc <- function(truth, prob) {
  if (length(unique(truth)) < 2 || all(is.na(prob))) return(NA_real_)
  suppressMessages(as.numeric(pROC::auc(
    pROC::roc(truth, prob, levels = c("Control", "CD"), direction = "<", quiet = TRUE))))
}

idx_ctrl <- which(dd$group == "Control")
idx_cd   <- which(dd$group == "CD")

message("\nRunning ", N_BOOT, " bootstrap replicates ...")
set.seed(2024)
boot <- vector("list", N_BOOT)
pb <- utils::txtProgressBar(min = 0, max = N_BOOT, style = 3)
for (b in seq_len(N_BOOT)) {
  # stratified resampling preserves the CD:Control ratio in every replicate
  ib <- c(sample(idx_ctrl, length(idx_ctrl), replace = TRUE),
          sample(idx_cd,   length(idx_cd),   replace = TRUE))
  oob <- setdiff(seq_len(nrow(dd)), unique(ib))

  train <- dd[ib, , drop = FALSE]
  if (length(unique(train$group)) < 2) { boot[[b]] <- NULL; next }
  m <- try(suppressWarnings(stats::glm(fml, data = train, family = stats::binomial())),
           silent = TRUE)
  if (inherits(m, "try-error")) { boot[[b]] <- NULL; next }

  auc_boot_in  <- safe_auc(train$group, stats::fitted(m))                      # apparent in resample
  auc_orig     <- safe_auc(dd$group, stats::predict(m, newdata = dd, type = "response"))
  auc_oob      <- if (length(oob) >= 4) {
    te <- dd[oob, , drop = FALSE]
    safe_auc(te$group, stats::predict(m, newdata = te, type = "response"))
  } else NA_real_

  boot[[b]] <- data.frame(replicate = b, n_oob = length(oob),
                          auc_oob = auc_oob, auc_boot_apparent = auc_boot_in,
                          auc_full_data = auc_orig)
  if (b %% 50 == 0) utils::setTxtProgressBar(pb, b)
}
close(pb)

reps <- do.call(rbind, Filter(Negate(is.null), boot))
data.table::fwrite(reps, file.path(RES_DIR, "bootstrap_auc_replicates.csv"))
message("Valid replicates: ", nrow(reps), "/", N_BOOT)

## ---------------------------------------------------------------------------
## 4. Summaries: percentile CI and optimism correction
## ---------------------------------------------------------------------------
alpha  <- (1 - CI_LEVEL) / 2
oob    <- reps$auc_oob[is.finite(reps$auc_oob)]
mean_oob <- mean(oob)
ci_oob   <- as.numeric(stats::quantile(oob, c(alpha, 1 - alpha), na.rm = TRUE))
se_oob   <- stats::sd(oob) / sqrt(length(oob))

optimism <- mean(reps$auc_boot_apparent - reps$auc_full_data, na.rm = TRUE)
auc_corrected <- auc_app - optimism
w632 <- 0.368 * auc_app + 0.632 * mean_oob    # Efron-Tibshirani .632 estimator

summary_tab <- data.frame(
  Metric = c("Apparent AUC (full data)",
             "Bootstrap out-of-bag AUC (mean)",
             "Bootstrap OOB 95% CI lower",
             "Bootstrap OOB 95% CI upper",
             "Bootstrap SE of the mean",
             "Optimism (Harrell)",
             "Optimism-corrected AUC",
             ".632 estimator",
             "Replicates",
             "Reference AUC (manuscript)",
             "Reference 95% CI (manuscript)"),
  Value = c(round(auc_app, 4), round(mean_oob, 4), round(ci_oob[1], 4), round(ci_oob[2], 4),
            round(se_oob, 5), round(optimism, 4), round(auc_corrected, 4), round(w632, 4),
            nrow(reps), REF_MEAN_AUC, paste0(REF_CI[1], " - ", REF_CI[2])),
  stringsAsFactors = FALSE)
data.table::fwrite(summary_tab, file.path(RES_DIR, "bootstrap_auc_summary.csv"))
print(summary_tab)

message(sprintf("\nBootstrap-validated AUC = %.3f (95%% CI %.3f - %.3f), %d resamples",
                mean_oob, ci_oob[1], ci_oob[2], nrow(reps)))
message(sprintf("Manuscript reference    = %.3f (95%% CI %.3f - %.3f), 2000 resamples",
                REF_MEAN_AUC, REF_CI[1], REF_CI[2]))

## ---------------------------------------------------------------------------
## 5. Optional: external cohort bootstrap (GSE75214)
## ---------------------------------------------------------------------------
ci_rows <- list(data.frame(Cohort = "GSE16879 (discovery, OOB bootstrap)",
                           N = nrow(dd), AUC = mean_oob,
                           CI_low = ci_oob[1], CI_high = ci_oob[2],
                           N_resamples = nrow(reps), stringsAsFactors = FALSE))

val_expr <- file.path(DATA_DIR, "GSE75214_expr_normalized.txt")
if (file.exists(val_expr)) {
  dv <- try(build_df("GSE75214_expr_normalized.txt", "GSE75214_phenotype.csv", GENES),
            silent = TRUE)
  if (!inherits(dv, "try-error")) {
    pv <- stats::predict(full, newdata = dv, type = "response")
    roc_v <- pROC::roc(dv$group, pv, levels = c("Control", "CD"), direction = "<", quiet = TRUE)
    set.seed(2024)
    ci_v <- pROC::ci.auc(roc_v, method = "bootstrap", boot.n = N_BOOT, progress = "none")
    message(sprintf("Validation cohort AUC   = %.3f (95%% CI %.3f - %.3f), stratified bootstrap",
                    as.numeric(pROC::auc(roc_v)), ci_v[1], ci_v[3]))
    ci_rows[[2]] <- data.frame(Cohort = "GSE75214 (external, stratified bootstrap CI)",
                               N = nrow(dv), AUC = as.numeric(pROC::auc(roc_v)),
                               CI_low = ci_v[1], CI_high = ci_v[3],
                               N_resamples = N_BOOT, stringsAsFactors = FALSE)
  }
}
ci_tab <- do.call(rbind, ci_rows)
ci_tab[, 3:5] <- round(ci_tab[, 3:5], 4)
data.table::fwrite(ci_tab, file.path(RES_DIR, "bootstrap_ci_table.csv"))

## ---------------------------------------------------------------------------
## 6. Figures
## ---------------------------------------------------------------------------
df_plot <- data.frame(AUC = oob)
p_dist <- ggplot(df_plot, aes(x = AUC)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45,
                 fill = "#2E86AB", colour = "white", alpha = 0.85) +
  geom_density(colour = "grey20", linewidth = 0.6) +
  geom_vline(xintercept = mean_oob, colour = "#E94F37", linewidth = 0.9) +
  geom_vline(xintercept = ci_oob, colour = "#E94F37", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = mean_oob, y = Inf, vjust = 1.6, hjust = -0.05, size = 3.6,
           colour = "#E94F37",
           label = sprintf("mean AUC = %.3f\n95%% CI %.3f - %.3f",
                           mean_oob, ci_oob[1], ci_oob[2])) +
  labs(title = "Bootstrap validation of the CFI/PAQR5/KCNE3 signature",
       subtitle = sprintf("%d stratified bootstrap resamples, out-of-bag evaluation (GSE16879)",
                          nrow(reps)),
       x = "Out-of-bag AUC", y = "Density") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Bootstrap_AUC_distribution.png"), p_dist,
       width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(FIG_DIR, "Bootstrap_AUC_distribution.pdf"), p_dist, width = 7.5, height = 5.5)

p_forest <- ggplot(ci_tab, aes(x = AUC, y = Cohort)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y",
                width = 0.12, linewidth = 0.7, colour = "#2E86AB") +
  geom_point(size = 3.2, colour = "#2E86AB") +
  geom_text(aes(label = sprintf("%.3f (%.3f-%.3f)", AUC, CI_low, CI_high)),
            vjust = -1.1, size = 3.3) +
  coord_cartesian(xlim = c(0.4, 1.02)) +
  labs(title = "Signature discrimination with bootstrap confidence intervals",
       x = "AUC", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Bootstrap_AUC_forest.png"), p_forest,
       width = 8.5, height = 3.6, dpi = 300)
ggsave(file.path(FIG_DIR, "Bootstrap_AUC_forest.pdf"), p_forest, width = 8.5, height = 3.6)

message("\nDone. Tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)

## ---------------------------------------------------------------------------
## 7. Session information
## ---------------------------------------------------------------------------
sessionInfo()
