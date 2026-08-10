###############################################################################
# 08c_individual_gene_AUC_corrected.R
#
# Purpose: Calculate AUC and 95% CI for individual genes with CORRECT direction
#          AUC should be > 0.5 for diagnostic performance
###############################################################################

suppressPackageStartupMessages({
  library(pROC)
  library(data.table)
})

set.seed(2024)

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

GENES <- c("CFI", "PAQR5", "KCNE3")
N_BOOT <- 2000

## ---------------------------------------------------------------------------
## Data loading
## ---------------------------------------------------------------------------
read_expr_matrix <- function(path) {
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

e  <- read_expr_matrix(file.path(DATA_DIR, "GSE16879_expr_normalized.txt"))
ph <- read_phenotype(file.path(DATA_DIR, "GSE16879_phenotype.csv"))
common <- intersect(colnames(e), ph$sample)
e  <- e[, common, drop = FALSE]; ph <- ph[match(common, ph$sample), , drop = FALSE]

group <- factor(ph$group, levels = c("Control", "CD"))

message("Discovery cohort: ", length(group), " samples (CD = ", sum(group == "CD"),
        ", Control = ", sum(group == "Control"), ")")

## ---------------------------------------------------------------------------
## Calculate AUC with CORRECTED direction (AUC > 0.5)
## ---------------------------------------------------------------------------
results <- data.frame(
  Gene = character(),
  AUC_raw = numeric(),
  AUC_corrected = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  Direction = character(),
  stringsAsFactors = FALSE
)

for (gene in GENES) {
  expr <- as.numeric(e[gene, ])
  
  # Determine direction based on group medians
  median_cd <- median(expr[group == "CD"])
  median_ctrl <- median(expr[group == "Control"])
  direction <- ifelse(median_cd > median_ctrl, ">", "<")
  
  # Calculate ROC
  roc_obj <- pROC::roc(group, expr, levels = c("Control", "CD"), 
                       direction = direction, quiet = TRUE)
  auc_raw <- as.numeric(pROC::auc(roc_obj))
  
  # Correct AUC to be > 0.5
  auc_corrected <- ifelse(auc_raw < 0.5, 1 - auc_raw, auc_raw)
  
  # Bootstrap CI
  set.seed(2024)
  boot_aucs <- numeric(N_BOOT)
  idx_ctrl <- which(group == "Control")
  idx_cd <- which(group == "CD")
  
  for (b in 1:N_BOOT) {
    ib <- c(sample(idx_ctrl, length(idx_ctrl), replace = TRUE),
            sample(idx_cd, length(idx_cd), replace = TRUE))
    boot_group <- group[ib]
    boot_expr <- expr[ib]
    
    if (length(unique(boot_group)) < 2) next
    
    median_cd_b <- median(boot_expr[boot_group == "CD"])
    median_ctrl_b <- median(boot_expr[boot_group == "Control"])
    dir_b <- ifelse(median_cd_b > median_ctrl_b, ">", "<")
    
    roc_boot <- try(pROC::roc(boot_group, boot_expr, levels = c("Control", "CD"),
                              direction = dir_b, quiet = TRUE), silent = TRUE)
    if (!inherits(roc_boot, "try-error")) {
      auc_b <- as.numeric(pROC::auc(roc_boot))
      boot_aucs[b] <- ifelse(auc_b < 0.5, 1 - auc_b, auc_b)
    }
  }
  
  ci_lower <- quantile(boot_aucs, 0.025, na.rm = TRUE)
  ci_upper <- quantile(boot_aucs, 0.975, na.rm = TRUE)
  
  results <- rbind(results, data.frame(
    Gene = gene,
    AUC_raw = round(auc_raw, 3),
    AUC_corrected = round(auc_corrected, 2),
    CI_lower = round(ci_lower, 2),
    CI_upper = round(ci_upper, 2),
    Direction = direction,
    stringsAsFactors = FALSE
  ))
  
  message(sprintf("%s: raw AUC = %.3f, corrected = %.2f (95%% CI %.2f - %.2f), direction = %s", 
                  gene, auc_raw, auc_corrected, ci_lower, ci_upper, direction))
}

message("\n=== CORRECTED Individual Gene AUC (GSE16879 Discovery) ===")
print(results[, c("Gene", "AUC_corrected", "CI_lower", "CI_upper")])

# Save
RES_DIR <- file.path(ROOT, "results", "08_bootstrap")
data.table::fwrite(results, file.path(RES_DIR, "individual_gene_auc_corrected.csv"))
message("\nSaved to: ", file.path(RES_DIR, "individual_gene_auc_corrected.csv"))
