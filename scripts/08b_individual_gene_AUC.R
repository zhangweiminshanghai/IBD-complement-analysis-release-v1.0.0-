###############################################################################
# 08b_individual_gene_AUC.R
#
# Purpose: Calculate AUC and 95% CI for individual genes (CFI, PAQR5, KCNE3)
#          in the discovery cohort (GSE16879) for Figure S6 Panel A
#
# Output: AUC with 95% CI for each gene
###############################################################################

suppressPackageStartupMessages({
  library(pROC)
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

GENES <- c("CFI", "PAQR5", "KCNE3")
N_BOOT <- 2000

## ---------------------------------------------------------------------------
## 1. Data loading
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

e  <- read_expr_matrix(file.path(DATA_DIR, "GSE16879_expr_normalized.txt"))
ph <- read_phenotype(file.path(DATA_DIR, "GSE16879_phenotype.csv"))
common <- intersect(colnames(e), ph$sample)
e  <- e[, common, drop = FALSE]; ph <- ph[match(common, ph$sample), , drop = FALSE]

message("Discovery cohort: ", nrow(ph), " samples (CD = ", sum(ph$group == "CD"),
        ", Control = ", sum(ph$group == "Control"), ")")

## ---------------------------------------------------------------------------
## 2. Calculate AUC and bootstrap CI for each gene
## ---------------------------------------------------------------------------
results <- data.frame(
  Gene = character(),
  AUC = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  N_CD = integer(),
  N_Control = integer(),
  stringsAsFactors = FALSE
)

for (gene in GENES) {
  if (!(gene %in% rownames(e))) {
    message("WARNING: ", gene, " not found in expression matrix")
    next
  }
  
  expr <- as.numeric(e[gene, ])
  group <- factor(ph$group, levels = c("Control", "CD"))
  
  # Calculate ROC and AUC
  roc_obj <- pROC::roc(group, expr, levels = c("Control", "CD"), 
                       direction = ifelse(median(expr[group == "CD"]) > median(expr[group == "Control"]), 
                                          ">", "<"), 
                       quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
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
    
    roc_boot <- try(pROC::roc(boot_group, boot_expr, levels = c("Control", "CD"),
                              direction = ifelse(median(boot_expr[boot_group == "CD"]) > 
                                                 median(boot_expr[boot_group == "Control"]), 
                                                ">", "<"),
                              quiet = TRUE), silent = TRUE)
    if (!inherits(roc_boot, "try-error")) {
      boot_aucs[b] <- as.numeric(pROC::auc(roc_boot))
    }
  }
  
  ci_lower <- quantile(boot_aucs, 0.025, na.rm = TRUE)
  ci_upper <- quantile(boot_aucs, 0.975, na.rm = TRUE)
  
  results <- rbind(results, data.frame(
    Gene = gene,
    AUC = round(auc_val, 2),
    CI_lower = round(ci_lower, 2),
    CI_upper = round(ci_upper, 2),
    N_CD = sum(group == "CD"),
    N_Control = sum(group == "Control"),
    stringsAsFactors = FALSE
  ))
  
  message(sprintf("%s: AUC = %.2f (95%% CI %.2f - %.2f)", 
                  gene, auc_val, ci_lower, ci_upper))
}

## ---------------------------------------------------------------------------
## 3. Output results
## ---------------------------------------------------------------------------
message("\n=== Individual Gene AUC Results (GSE16879 Discovery) ===")
print(results)

# Save results
RES_DIR <- file.path(ROOT, "results", "08_bootstrap")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(results, file.path(RES_DIR, "individual_gene_auc.csv"))

message("\nResults saved to: ", file.path(RES_DIR, "individual_gene_auc.csv"))
