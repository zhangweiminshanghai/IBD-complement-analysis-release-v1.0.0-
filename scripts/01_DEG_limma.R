###############################################################################
# 01_DEG_limma.R
#
# Purpose      : Differential expression analysis (limma) of Crohn's disease (CD)
#                vs non-IBD control intestinal mucosa in two independent cohorts:
#                  * GSE16879  (discovery,  73 CD / 12 controls, hgu133plus2)
#                  * GSE75214  (validation, 75 CD / 22 controls, hugene10sttranscriptcluster)
#                Produces per-cohort DEG tables, the shared DEG list, and volcano plots.
#
# Inputs       : data/processed/GSE16879_expr_normalized.txt
#                    tab-delimited; column 1 = "Gene Symbol" (or probe ID), remaining
#                    columns = GSM sample IDs. GSE16879 is stored on the LINEAR scale
#                    (MAS5-like) and is log2-transformed automatically.
#                data/processed/GSE16879_phenotype.csv
#                    CSV with at least: sample,group   (group in {CD, Control})
#                    optional columns: tissue, timepoint, response, batch
#                    NOTE ON GSE16879: the series profiles mucosa BEFORE and AFTER
#                    infliximab. The published discovery contrast uses baseline
#                    (pre-treatment) CD mucosa vs non-IBD controls. When the
#                    phenotype table carries a timepoint/treatment column,
#                    RESTRICT_TO_BASELINE = TRUE keeps week-0 CD samples only
#                    (controls are always retained). Mixing post-treatment
#                    responders into the CD arm dilutes the contrast heavily.
#                data/processed/GSE75214_expr_normalized.txt   (already log2)
#                data/processed/GSE75214_phenotype.csv
#
# Outputs      : results/01_DEG/GSE16879_DEG_all.csv
#                results/01_DEG/GSE16879_DEG_significant.csv
#                results/01_DEG/GSE75214_DEG_all.csv
#                results/01_DEG/GSE75214_DEG_significant.csv
#                results/01_DEG/Common_DEGs.csv
#                results/01_DEG/Common_DEGs_for_STRING.txt
#                results/01_DEG/DEG_summary.csv
#                results/figures/VolcanoPlot_GSE16879.{png,pdf}
#                results/figures/VolcanoPlot_GSE75214.{png,pdf}
#                results/figures/Combined_VolcanoPlots_CD_vs_Control.{png,pdf}
#
# Figure/Table : Figure 1A-B (volcano plots, discovery + validation)
#                Supplementary Table S2 (full DEG tables)
#
# Key params   : DEG thresholds |log2FC| > 1 AND adj.P.Val < 0.05 (BH)
#                limma design: ~ 0 + group ; contrast = CD - Control
#                set.seed(2024)
#
# Runtime      : ~2-4 min on a laptop (both cohorts, incl. plotting)
#
# Packages     : limma (>= 3.54), ggplot2 (>= 3.4), ggrepel (>= 0.9.3),
#                patchwork (>= 1.1.2), data.table (>= 1.14)
#                Optional (only if the input matrix is still probe-level):
#                  hgu133plus2.db (>= 3.13), hugene10sttranscriptcluster.db (>= 8.8),
#                  AnnotationDbi (>= 1.60)
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(data.table)
})

set.seed(2024)

## ---------------------------------------------------------------------------
## 0. Paths -- resolved relative to the repository root so the script runs both
##    from the repo root and from inside scripts/
## ---------------------------------------------------------------------------
find_repo_root <- function() {
  cand <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:4) {
    if (dir.exists(file.path(cand, "data")) && dir.exists(file.path(cand, "scripts"))) {
      return(cand)
    }
    cand <- normalizePath(file.path(cand, ".."), winslash = "/", mustWork = FALSE)
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

ROOT     <- find_repo_root()
DATA_DIR <- file.path(ROOT, "data", "processed")
RES_DIR  <- file.path(ROOT, "results", "01_DEG")
FIG_DIR  <- file.path(ROOT, "results", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LOGFC_CUT <- 1
PADJ_CUT  <- 0.05
# Keep pre-treatment CD samples when the phenotype table documents a
# timepoint / treatment column (see the note in the header).
RESTRICT_TO_BASELINE <- TRUE

message("Repository root: ", ROOT)

## ---------------------------------------------------------------------------
## 1. Helpers
## ---------------------------------------------------------------------------

#' Read a processed expression matrix (feature x sample) from a tab-delimited file.
#' Duplicated feature names are collapsed by keeping the row with the highest mean.
read_expr_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Missing expression matrix: ", path,
         "\n  -> place the processed matrix in data/processed/ (see header).")
  }
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  feat <- as.character(dt[[1]])
  mat  <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feat

  mat <- mat[!is.na(feat) & feat != "" & feat != "---", , drop = FALSE]
  if (anyDuplicated(rownames(mat))) {
    ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
    mat <- mat[ord, , drop = FALSE]
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }
  mat[order(rownames(mat)), , drop = FALSE]
}

#' log2-transform if the matrix is clearly on the linear scale (GSE16879 case).
maybe_log2 <- function(mat, label = "") {
  qx <- as.numeric(stats::quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  needs_log <- (qx[6] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[5] > 2)
  if (needs_log) {
    message("  [", label, "] linear scale detected (max = ", round(qx[6], 1),
            ") -> applying log2(x + 1)")
    mat[mat < 0] <- 0
    mat <- log2(mat + 1)
  } else {
    message("  [", label, "] already log2-scaled (max = ", round(qx[6], 2), ")")
  }
  mat
}

#' Probe -> gene symbol collapse. Only used when the matrix is still probe-level.
#' GSE16879 : hgu133plus2SYMBOL ; GSE75214 : hugene10sttranscriptclusterSYMBOL
probe2symbol <- function(mat, platform = c("hgu133plus2", "hugene10sttranscriptcluster")) {
  platform <- match.arg(platform)
  looks_like_symbol <- mean(grepl("^[A-Za-z][A-Za-z0-9._-]*$", rownames(mat))) > 0.8 &&
    mean(grepl("_at$|^\\d+$", rownames(mat))) < 0.2
  if (looks_like_symbol) {
    message("  rownames already look like gene symbols - skipping probe annotation")
    return(mat)
  }
  pkg <- paste0(platform, ".db")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Probe-level matrix supplied but annotation package '", pkg, "' is not installed.\n",
         "  install: BiocManager::install('", pkg, "')")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  map    <- get(paste0(platform, "SYMBOL"))
  ids    <- AnnotationDbi::mappedkeys(map)
  p2s    <- unlist(as.list(map[ids]))
  keep   <- intersect(rownames(mat), names(p2s))
  message("  probe annotation: ", length(keep), "/", nrow(mat), " probes mapped")
  mat    <- mat[keep, , drop = FALSE]
  symbol <- p2s[keep]
  # collapse duplicated symbols using the maximum-mean probe
  ord    <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat    <- mat[ord, , drop = FALSE]
  symbol <- symbol[ord]
  keep2  <- !duplicated(symbol)
  mat    <- mat[keep2, , drop = FALSE]
  rownames(mat) <- symbol[keep2]
  mat[order(rownames(mat)), , drop = FALSE]
}

#' Read a phenotype table and align it to the expression matrix columns.
read_phenotype <- function(path, expr) {
  if (!file.exists(path)) {
    stop("Missing phenotype file: ", path,
         "\n  -> expected CSV with columns 'sample' and 'group' (CD / Control).")
  }
  ph <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  cn <- tolower(colnames(ph))
  scol <- colnames(ph)[which(cn %in% c("sample", "gsm", "geo_accession", "sample_id"))[1]]
  gcol <- colnames(ph)[which(cn %in% c("group", "class", "condition", "diagnosis", "disease"))[1]]
  if (is.na(scol) || is.na(gcol)) {
    stop("Phenotype file must contain a sample column and a group column: ", path)
  }
  ph$sample <- as.character(ph[[scol]])
  ph$group  <- as.character(ph[[gcol]])
  # harmonise labels
  ph$group[grepl("^(cd|crohn)", ph$group, ignore.case = TRUE)]                 <- "CD"
  ph$group[grepl("^(control|normal|non.?ibd|healthy)", ph$group, ignore.case = TRUE)] <- "Control"
  ph <- ph[ph$group %in% c("CD", "Control"), , drop = FALSE]   # drop UC / other

  # optional baseline restriction (GSE16879 before/after infliximab design)
  tcol <- colnames(ph)[which(tolower(colnames(ph)) %in%
                               c("timepoint", "time", "treatment", "visit", "week"))[1]]
  if (RESTRICT_TO_BASELINE && !is.na(tcol)) {
    baseline <- grepl("before|baseline|pre|w0|week ?0|untreated", ph[[tcol]], ignore.case = TRUE)
    drop <- ph$group == "CD" & !baseline
    if (any(drop)) {
      message("  restricting CD arm to baseline samples: dropping ", sum(drop),
              " post-treatment sample(s) using column '", tcol, "'")
      ph <- ph[!drop, , drop = FALSE]
    }
  }

  common <- intersect(colnames(expr), ph$sample)
  if (length(common) < 10) {
    stop("Fewer than 10 samples shared between expression matrix and phenotype table.")
  }
  ph <- ph[match(common, ph$sample), , drop = FALSE]
  list(expr = expr[, common, drop = FALSE], pheno = ph)
}

#' limma differential expression: CD vs Control.
run_limma <- function(expr, group) {
  group <- factor(group, levels = c("Control", "CD"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit  <- limma::lmFit(as.matrix(expr), design)
  cont <- limma::makeContrasts(CDvsControl = CD - Control, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  res  <- limma::topTable(fit2, coef = "CDvsControl", number = Inf, sort.by = "P")
  res$Gene <- rownames(res)
  res$Change <- "NoChange"
  res$Change[res$logFC >  LOGFC_CUT & res$adj.P.Val < PADJ_CUT] <- "Up"
  res$Change[res$logFC < -LOGFC_CUT & res$adj.P.Val < PADJ_CUT] <- "Down"
  res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Change")]
}

#' Publication-style volcano plot with the top-10 genes labelled.
volcano_plot <- function(res, title, subtitle = NULL, n_label = 10) {
  df <- res
  df$Change <- factor(df$Change, levels = c("Down", "NoChange", "Up"))
  top <- head(df[df$Change != "NoChange", ][order(df$P.Value[df$Change != "NoChange"]), ], n_label)
  # always show the study's core genes when they pass the filter
  core <- df[df$Gene %in% c("CFI", "PAQR5", "KCNE3", "C2", "CFB") & df$Change != "NoChange", ]
  lab  <- unique(rbind(top, core))

  ggplot(df, aes(x = logFC, y = -log10(P.Value), colour = Change)) +
    geom_point(alpha = 0.55, size = 1.2) +
    scale_colour_manual(values = c(Down = "#2E86AB", NoChange = "grey78", Up = "#E94F37"),
                        labels = c(
                          Down     = paste0("Down (", sum(df$Change == "Down"), ")"),
                          NoChange = "NS",
                          Up       = paste0("Up (", sum(df$Change == "Up"), ")"))) +
    geom_vline(xintercept = c(-LOGFC_CUT, LOGFC_CUT), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(max(df$P.Value[df$adj.P.Val < PADJ_CUT], 1e-300)),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = lab, aes(label = Gene), size = 3,
                             max.overlaps = 30, box.padding = 0.35,
                             segment.size = 0.25, show.legend = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2] * " fold change (CD vs Control)"),
         y = expression(-log[10] * " (P value)"), colour = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
          plot.title = element_text(face = "bold"))
}

## ---------------------------------------------------------------------------
## 2. Cohort configuration
## ---------------------------------------------------------------------------
cohorts <- list(
  GSE16879 = list(
    expr_file  = file.path(DATA_DIR, "GSE16879_expr_normalized.txt"),
    pheno_file = file.path(DATA_DIR, "GSE16879_phenotype.csv"),
    platform   = "hgu133plus2",
    label      = "GSE16879 (discovery: 73 CD vs 12 control)"
  ),
  GSE75214 = list(
    expr_file  = file.path(DATA_DIR, "GSE75214_expr_normalized.txt"),
    pheno_file = file.path(DATA_DIR, "GSE75214_phenotype.csv"),
    platform   = "hugene10sttranscriptcluster",
    label      = "GSE75214 (validation: 75 CD vs 22 control)"
  )
)

## ---------------------------------------------------------------------------
## 3. Run limma per cohort
## ---------------------------------------------------------------------------
deg_list  <- list()
plot_list <- list()
summary_rows <- list()

for (nm in names(cohorts)) {
  cfg <- cohorts[[nm]]
  message("\n=== ", nm, " ===")

  if (!file.exists(cfg$expr_file) || !file.exists(cfg$pheno_file)) {
    warning("Skipping ", nm, ": missing ",
            paste(basename(c(cfg$expr_file, cfg$pheno_file)[
              !file.exists(c(cfg$expr_file, cfg$pheno_file))]), collapse = " and "))
    next
  }

  expr <- read_expr_matrix(cfg$expr_file)
  expr <- maybe_log2(expr, nm)
  expr <- probe2symbol(expr, cfg$platform)

  al    <- read_phenotype(cfg$pheno_file, expr)
  expr  <- al$expr
  pheno <- al$pheno
  message("  matrix: ", nrow(expr), " genes x ", ncol(expr), " samples  (CD = ",
          sum(pheno$group == "CD"), ", Control = ", sum(pheno$group == "Control"), ")")

  # drop uninformative rows (no variance / all NA)
  keep <- apply(expr, 1, function(x) sum(!is.na(x)) > 3 && stats::sd(x, na.rm = TRUE) > 0)
  expr <- expr[keep, , drop = FALSE]

  res <- run_limma(expr, pheno$group)
  deg_list[[nm]] <- res

  sig <- res[res$Change != "NoChange", ]
  data.table::fwrite(res, file.path(RES_DIR, paste0(nm, "_DEG_all.csv")))
  data.table::fwrite(sig, file.path(RES_DIR, paste0(nm, "_DEG_significant.csv")))
  message("  DEGs: up = ", sum(res$Change == "Up"),
          " | down = ", sum(res$Change == "Down"),
          " | total = ", nrow(sig))

  summary_rows[[nm]] <- data.frame(
    cohort = nm, n_genes = nrow(res),
    n_CD = sum(pheno$group == "CD"), n_Control = sum(pheno$group == "Control"),
    n_up = sum(res$Change == "Up"), n_down = sum(res$Change == "Down"),
    n_DEG = nrow(sig),
    CFI_logFC   = res$logFC[match("CFI", res$Gene)],
    CFI_adjP    = res$adj.P.Val[match("CFI", res$Gene)],
    C2_logFC    = res$logFC[match("C2", res$Gene)],
    C2_adjP     = res$adj.P.Val[match("C2", res$Gene)],
    PAQR5_logFC = res$logFC[match("PAQR5", res$Gene)],
    KCNE3_logFC = res$logFC[match("KCNE3", res$Gene)],
    stringsAsFactors = FALSE
  )

  p <- volcano_plot(res, title = nm, subtitle = cfg$label)
  plot_list[[nm]] <- p
  ggsave(file.path(FIG_DIR, paste0("VolcanoPlot_", nm, ".png")), p,
         width = 7, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0("VolcanoPlot_", nm, ".pdf")), p,
         width = 7, height = 6)
}

## ---------------------------------------------------------------------------
## 4. Combined volcano figure (Figure 1A-B)
## ---------------------------------------------------------------------------
if (length(plot_list) == 2) {
  combined <- (plot_list[[1]] | plot_list[[2]]) +
    patchwork::plot_annotation(
      title = "Differential expression in Crohn's disease mucosa",
      subtitle = sprintf("limma, |log2FC| > %g and BH-adjusted P < %g", LOGFC_CUT, PADJ_CUT),
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(FIG_DIR, "Combined_VolcanoPlots_CD_vs_Control.png"), combined,
         width = 13, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, "Combined_VolcanoPlots_CD_vs_Control.pdf"), combined,
         width = 13, height = 6)
}

## ---------------------------------------------------------------------------
## 5. Shared DEGs across discovery and validation
## ---------------------------------------------------------------------------
if (length(deg_list) == 2) {
  d1 <- deg_list[["GSE16879"]]; d2 <- deg_list[["GSE75214"]]
  s1 <- d1[d1$Change != "NoChange", ]; s2 <- d2[d2$Change != "NoChange", ]
  common <- intersect(s1$Gene, s2$Gene)
  cm <- data.frame(
    Gene            = common,
    logFC_GSE16879  = s1$logFC[match(common, s1$Gene)],
    adjP_GSE16879   = s1$adj.P.Val[match(common, s1$Gene)],
    logFC_GSE75214  = s2$logFC[match(common, s2$Gene)],
    adjP_GSE75214   = s2$adj.P.Val[match(common, s2$Gene)],
    stringsAsFactors = FALSE
  )
  cm$Concordant <- sign(cm$logFC_GSE16879) == sign(cm$logFC_GSE75214)
  cm <- cm[order(cm$adjP_GSE16879), ]
  data.table::fwrite(cm, file.path(RES_DIR, "Common_DEGs.csv"))
  writeLines(cm$Gene[cm$Concordant], file.path(RES_DIR, "Common_DEGs_for_STRING.txt"))
  message("\nShared DEGs: ", nrow(cm), " (concordant direction: ", sum(cm$Concordant), ")")
}

data.table::fwrite(do.call(rbind, summary_rows), file.path(RES_DIR, "DEG_summary.csv"))

message("\nDone. Tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)

## ---------------------------------------------------------------------------
## 6. Session information (reproducibility)
## ---------------------------------------------------------------------------
sessionInfo()
