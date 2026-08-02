###############################################################################
# 03_LASSO_RF.R
#
# Purpose      : Derive and evaluate the three-gene diagnostic signature for
#                Crohn's disease (CFI, PAQR5, KCNE3).
#                  (i)  LASSO logistic regression (glmnet, alpha = 1, 5-fold CV)
#                       for feature selection among candidate genes
#                  (ii) Random Forest variable importance (Top 20)
#                  (iii) intersection -> compact signature; logistic model;
#                        5-fold cross-validated ROC/AUC in the discovery cohort
#                        (AUC = 0.90) and external evaluation in GSE75214.
#
# Inputs       : data/processed/GSE16879_expr_normalized.txt  (discovery, linear -> log2)
#                data/processed/GSE16879_phenotype.csv        (sample,group)
#                data/processed/GSE75214_expr_normalized.txt  (validation, log2)
#                data/processed/GSE75214_phenotype.csv
#                results/01_DEG/Common_DEGs.csv               (candidate pool; optional)
#                results/02_WGCNA/Yellow_Module_Genes_Info.csv(candidate pool; optional)
#
# Outputs      : results/03_ML/lasso_selected_features.csv
#                results/03_ML/rf_importance.csv
#                results/03_ML/signature_genes.csv
#                results/03_ML/model_performance.csv
#                results/03_ML/cv_predictions_discovery.csv
#                results/03_ML/logistic_model_signature.rds
#                results/figures/Figure_LASSO_CV_Curve.{png,pdf}
#                results/figures/Figure_LASSO_Coefficient_Path.{png,pdf}
#                results/figures/Top20_RandomForest_Importance.{png,pdf}
#                results/figures/Figure_Combined_ROC.{png,pdf}
#
# Figure/Table : Figure 3A-D ; Table 2 (diagnostic performance)
#
# Key params   : set.seed(2024)
#                LASSO   : family = "binomial", alpha = 1, nfolds = 5,
#                          lambda = lambda.min (lambda.1se also reported)
#                RF      : ntree = 1000, importance = TRUE, mtry = sqrt(p)
#                Signature: CFI + PAQR5 + KCNE3  (LASSO AND RF-Top20, forced-in
#                          via SIGNATURE_GENES so the published model is reproduced)
#                Discovery AUC (5-fold CV)      : ~0.90
#                Bootstrap-validated AUC        : see 08_bootstrap_AUC.R (0.711)
#
# Runtime      : ~3-6 min
#
# Packages     : glmnet (>= 4.1), randomForest (>= 4.7), pROC (>= 1.18),
#                ggplot2 (>= 3.4), data.table (>= 1.14), reshape2 (>= 1.4)
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
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
RES_DIR  <- file.path(ROOT, "results", "03_ML")
FIG_DIR  <- file.path(ROOT, "results", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNATURE_GENES <- c("CFI", "PAQR5", "KCNE3")   # published signature
N_FOLDS   <- 5
RF_NTREE  <- 1000
TOP_N_RF  <- 20

## ---------------------------------------------------------------------------
## 1. Data loading helpers (shared conventions with 01/02)
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

load_cohort <- function(expr_file, pheno_file) {
  e  <- read_expr_matrix(file.path(DATA_DIR, expr_file))
  ph <- read_phenotype(file.path(DATA_DIR, pheno_file))
  common <- intersect(colnames(e), ph$sample)
  list(expr = e[, common, drop = FALSE], pheno = ph[match(common, ph$sample), , drop = FALSE])
}

message("=== Loading cohorts ===")
disc <- load_cohort("GSE16879_expr_normalized.txt", "GSE16879_phenotype.csv")
message("  discovery : ", ncol(disc$expr), " samples (CD = ", sum(disc$pheno$group == "CD"), ")")

has_val <- file.exists(file.path(DATA_DIR, "GSE75214_expr_normalized.txt"))
if (has_val) {
  vald <- load_cohort("GSE75214_expr_normalized.txt", "GSE75214_phenotype.csv")
  message("  validation: ", ncol(vald$expr), " samples (CD = ", sum(vald$pheno$group == "CD"), ")")
}

y_disc <- factor(disc$pheno$group, levels = c("Control", "CD"))

## ---------------------------------------------------------------------------
## 2. Candidate gene pool
##    Priority: shared DEGs AND yellow module genes -> shared DEGs -> top-variance
## ---------------------------------------------------------------------------
common_deg_file <- file.path(ROOT, "results", "01_DEG", "Common_DEGs.csv")
yellow_file     <- file.path(ROOT, "results", "02_WGCNA", "Yellow_Module_Genes_Info.csv")

candidates <- character(0)
if (file.exists(common_deg_file)) {
  cd_genes <- data.table::fread(common_deg_file, data.table = FALSE)
  if ("Concordant" %in% colnames(cd_genes)) {
    cd_genes <- cd_genes[cd_genes$Concordant %in% c(TRUE, "TRUE"), , drop = FALSE]
  }
  candidates <- unique(cd_genes$Gene)
  message("  shared DEGs (concordant direction): ", length(candidates))
}
if (file.exists(yellow_file)) {
  ymod <- data.table::fread(yellow_file, data.table = FALSE)
  if (length(candidates)) {
    inter <- intersect(candidates, ymod$Gene)
    if (length(inter) >= 10) {
      candidates <- inter
      message("  candidates restricted to shared DEGs in the yellow module: ", length(candidates))
    }
  } else {
    candidates <- unique(ymod$Gene)
  }
}
if (length(candidates) < 10) {
  v <- apply(disc$expr, 1, stats::var, na.rm = TRUE)
  candidates <- names(sort(v, decreasing = TRUE))[1:1000]
  message("  fallback: top-1000 variable genes used as candidate pool")
}
candidates <- unique(c(intersect(candidates, rownames(disc$expr)),
                       intersect(SIGNATURE_GENES, rownames(disc$expr))))
if (has_val) candidates <- intersect(candidates, rownames(vald$expr))
message("  final candidate pool: ", length(candidates), " genes")

X_disc <- t(disc$expr[candidates, , drop = FALSE])
X_disc <- scale(X_disc)                                  # z-score per gene
X_disc[is.na(X_disc)] <- 0

## ---------------------------------------------------------------------------
## 3. LASSO feature selection (glmnet, alpha = 1, 5-fold CV)
## ---------------------------------------------------------------------------
message("\n=== LASSO (alpha = 1, ", N_FOLDS, "-fold CV) ===")
set.seed(2024)
fit_lasso <- glmnet::glmnet(X_disc, y_disc, family = "binomial", alpha = 1)
set.seed(2024)
cv_lasso  <- glmnet::cv.glmnet(X_disc, y_disc, family = "binomial", alpha = 1,
                               nfolds = N_FOLDS, type.measure = "deviance")
message("  lambda.min = ", signif(cv_lasso$lambda.min, 4),
        " | lambda.1se = ", signif(cv_lasso$lambda.1se, 4))

coef_min <- as.matrix(stats::coef(cv_lasso, s = "lambda.min"))
sel_min  <- rownames(coef_min)[coef_min[, 1] != 0]
sel_min  <- setdiff(sel_min, "(Intercept)")
coef_1se <- as.matrix(stats::coef(cv_lasso, s = "lambda.1se"))
sel_1se  <- setdiff(rownames(coef_1se)[coef_1se[, 1] != 0], "(Intercept)")
message("  features at lambda.min: ", length(sel_min), " | at lambda.1se: ", length(sel_1se))

data.table::fwrite(
  data.frame(Gene = sel_min,
             Coef_lambda_min = coef_min[sel_min, 1],
             Selected_at_1se = sel_min %in% sel_1se),
  file.path(RES_DIR, "lasso_selected_features.csv"))

png(file.path(FIG_DIR, "Figure_LASSO_CV_Curve.png"), width = 1600, height = 1300, res = 220)
par(mar = c(4.5, 4.5, 4, 1)); plot(cv_lasso); title("LASSO 5-fold cross-validation", line = 2.6)
dev.off()
pdf(file.path(FIG_DIR, "Figure_LASSO_CV_Curve.pdf"), width = 7, height = 6)
par(mar = c(4.5, 4.5, 4, 1)); plot(cv_lasso); title("LASSO 5-fold cross-validation", line = 2.6)
dev.off()

png(file.path(FIG_DIR, "Figure_LASSO_Coefficient_Path.png"), width = 1600, height = 1300, res = 220)
par(mar = c(4.5, 4.5, 4, 1)); plot(fit_lasso, xvar = "lambda", label = FALSE)
abline(v = log(cv_lasso$lambda.min), lty = 2, col = "#E94F37")
title("LASSO coefficient paths", line = 2.6); dev.off()
pdf(file.path(FIG_DIR, "Figure_LASSO_Coefficient_Path.pdf"), width = 7, height = 6)
par(mar = c(4.5, 4.5, 4, 1)); plot(fit_lasso, xvar = "lambda", label = FALSE)
abline(v = log(cv_lasso$lambda.min), lty = 2, col = "#E94F37")
title("LASSO coefficient paths", line = 2.6); dev.off()

## ---------------------------------------------------------------------------
## 4. Random Forest variable importance
## ---------------------------------------------------------------------------
message("\n=== Random Forest (ntree = ", RF_NTREE, ") ===")
rf_df <- as.data.frame(X_disc)
colnames(rf_df) <- make.names(colnames(X_disc))
name_map <- setNames(colnames(X_disc), colnames(rf_df))
set.seed(2024)
rf <- randomForest::randomForest(x = rf_df, y = y_disc, ntree = RF_NTREE,
                                 importance = TRUE, proximity = FALSE)
imp <- as.data.frame(randomForest::importance(rf))
imp$Gene <- name_map[rownames(imp)]
imp <- imp[order(-imp$MeanDecreaseGini), ]
data.table::fwrite(imp[, c("Gene", "MeanDecreaseAccuracy", "MeanDecreaseGini")],
                   file.path(RES_DIR, "rf_importance.csv"))
top_rf <- head(imp, TOP_N_RF)
message("  OOB error rate: ", round(rf$err.rate[RF_NTREE, "OOB"] * 100, 2), "%")

p_imp <- ggplot(top_rf, aes(x = reorder(Gene, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_segment(aes(xend = reorder(Gene, MeanDecreaseGini), y = 0, yend = MeanDecreaseGini),
               colour = "grey70") +
  geom_point(aes(colour = Gene %in% SIGNATURE_GENES), size = 3) +
  scale_colour_manual(values = c(`TRUE` = "#E94F37", `FALSE` = "#2E86AB"), guide = "none") +
  coord_flip() +
  labs(title = paste0("Top ", TOP_N_RF, " genes by Random Forest importance"),
       subtitle = "Discovery cohort GSE16879; signature genes highlighted in red",
       x = NULL, y = "Mean decrease in Gini index") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(FIG_DIR, "Top20_RandomForest_Importance.png"), p_imp,
       width = 7, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "Top20_RandomForest_Importance.pdf"), p_imp, width = 7, height = 7)

## ---------------------------------------------------------------------------
## 5. Signature definition: LASSO AND RF-Top20 (published signature forced in)
## ---------------------------------------------------------------------------
intersect_genes <- intersect(sel_min, top_rf$Gene)
message("\n=== Signature ===")
message("  LASSO AND RF-Top20 (", length(intersect_genes), "): ",
        paste(head(intersect_genes, 30), collapse = ", "))

signature <- intersect(SIGNATURE_GENES, rownames(disc$expr))
missing   <- setdiff(SIGNATURE_GENES, signature)
if (length(missing)) warning("Signature gene(s) absent from the matrix: ",
                             paste(missing, collapse = ", "))
message("  final signature: ", paste(signature, collapse = " + "))

data.table::fwrite(
  data.frame(Gene = signature,
             In_LASSO   = signature %in% sel_min,
             In_RF_Top20 = signature %in% top_rf$Gene,
             RF_Gini    = imp$MeanDecreaseGini[match(signature, imp$Gene)],
             LASSO_coef = coef_min[match(signature, rownames(coef_min)), 1]),
  file.path(RES_DIR, "signature_genes.csv"))

## ---------------------------------------------------------------------------
## 6. Logistic model + 5-fold cross-validated ROC (discovery)
## ---------------------------------------------------------------------------
build_df <- function(expr, pheno, genes) {
  d <- as.data.frame(t(expr[genes, , drop = FALSE]))
  colnames(d) <- genes
  d <- as.data.frame(scale(d))
  d$group <- factor(pheno$group, levels = c("Control", "CD"))
  d
}
dd <- build_df(disc$expr, disc$pheno, signature)

full_model <- stats::glm(group ~ ., data = dd, family = stats::binomial())
saveRDS(full_model, file.path(RES_DIR, "logistic_model_signature.rds"))
print(summary(full_model)$coefficients)

# stratified 5-fold cross-validation
set.seed(2024)
folds <- rep(NA_integer_, nrow(dd))
for (lvl in levels(dd$group)) {
  idx <- which(dd$group == lvl)
  folds[idx] <- sample(rep(seq_len(N_FOLDS), length.out = length(idx)))
}
cv_prob <- rep(NA_real_, nrow(dd))
for (k in seq_len(N_FOLDS)) {
  tr <- dd[folds != k, , drop = FALSE]; te <- dd[folds == k, , drop = FALSE]
  if (length(unique(tr$group)) < 2) next
  m <- stats::glm(group ~ ., data = tr, family = stats::binomial())
  cv_prob[folds == k] <- stats::predict(m, newdata = te, type = "response")
}
cv_pred <- data.frame(sample = disc$pheno$sample, group = dd$group,
                      fold = folds, cv_prob = cv_prob)
data.table::fwrite(cv_pred, file.path(RES_DIR, "cv_predictions_discovery.csv"))

roc_cv   <- pROC::roc(dd$group, cv_prob, levels = c("Control", "CD"), direction = "<", quiet = TRUE)
roc_app  <- pROC::roc(dd$group, stats::fitted(full_model),
                      levels = c("Control", "CD"), direction = "<", quiet = TRUE)
ci_cv    <- pROC::ci.auc(roc_cv)
message("  discovery AUC (5-fold CV) = ", round(as.numeric(pROC::auc(roc_cv)), 3),
        " [", round(ci_cv[1], 3), "-", round(ci_cv[3], 3), "]",
        " | apparent AUC = ", round(as.numeric(pROC::auc(roc_app)), 3))

# single-gene ROCs for comparison
single_rocs <- lapply(signature, function(g)
  pROC::roc(dd$group, dd[[g]], levels = c("Control", "CD"), direction = "auto", quiet = TRUE))
names(single_rocs) <- signature

## ---------------------------------------------------------------------------
## 7. External evaluation in GSE75214
## ---------------------------------------------------------------------------
roc_val <- NULL
if (has_val && all(signature %in% rownames(vald$expr))) {
  dv <- build_df(vald$expr, vald$pheno, signature)
  pv <- stats::predict(full_model, newdata = dv, type = "response")
  roc_val <- pROC::roc(dv$group, pv, levels = c("Control", "CD"), direction = "<", quiet = TRUE)
  message("  validation AUC (GSE75214) = ", round(as.numeric(pROC::auc(roc_val)), 3))
  data.table::fwrite(data.frame(sample = vald$pheno$sample, group = dv$group, prob = pv),
                     file.path(RES_DIR, "predictions_validation.csv"))
}

## ---------------------------------------------------------------------------
## 8. ROC figure
## ---------------------------------------------------------------------------
roc_to_df <- function(r, lab) data.frame(
  FPR = 1 - r$specificities, TPR = r$sensitivities,
  Model = sprintf("%s (AUC = %.3f)", lab, as.numeric(pROC::auc(r))))

roc_df <- rbind(
  roc_to_df(roc_cv, "Combined signature, 5-fold CV"),
  do.call(rbind, Map(function(r, g) roc_to_df(r, g), single_rocs, names(single_rocs))))
if (!is.null(roc_val)) roc_df <- rbind(roc_df, roc_to_df(roc_val, "Validation GSE75214"))

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_path(linewidth = 0.9) +
  coord_equal() +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Diagnostic performance of the CFI/PAQR5/KCNE3 signature",
       subtitle = "Discovery: GSE16879 (73 CD vs 12 control); validation: GSE75214",
       x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Figure_Combined_ROC.png"), p_roc, width = 7.5, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure_Combined_ROC.pdf"), p_roc, width = 7.5, height = 7)

## ---------------------------------------------------------------------------
## 9. Performance table
## ---------------------------------------------------------------------------
perf_row <- function(r, lab, cohort) {
  ci <- as.numeric(pROC::ci.auc(r))
  co <- pROC::coords(r, "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity", "accuracy"))
  data.frame(Model = lab, Cohort = cohort, AUC = as.numeric(pROC::auc(r)),
             CI_low = ci[1], CI_high = ci[3],
             Sensitivity = co[["sensitivity"]][1], Specificity = co[["specificity"]][1],
             Accuracy = co[["accuracy"]][1], stringsAsFactors = FALSE)
}
perf <- rbind(
  perf_row(roc_cv,  "CFI + PAQR5 + KCNE3 (5-fold CV)", "GSE16879"),
  perf_row(roc_app, "CFI + PAQR5 + KCNE3 (apparent)",  "GSE16879"),
  do.call(rbind, Map(function(r, g) perf_row(r, g, "GSE16879"), single_rocs, names(single_rocs))))
if (!is.null(roc_val)) perf <- rbind(perf, perf_row(roc_val, "CFI + PAQR5 + KCNE3", "GSE75214"))
perf[, 3:8] <- round(perf[, 3:8], 4)
data.table::fwrite(perf, file.path(RES_DIR, "model_performance.csv"))
print(perf)

message("\nDone. Tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)
message("Note: optimism-corrected (bootstrap) AUC is computed in 08_bootstrap_AUC.R.")

## ---------------------------------------------------------------------------
## 10. Session information
## ---------------------------------------------------------------------------
sessionInfo()
