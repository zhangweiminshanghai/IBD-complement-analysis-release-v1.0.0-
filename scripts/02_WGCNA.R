###############################################################################
# 02_WGCNA.R
#
# Purpose      : Weighted gene co-expression network analysis (WGCNA) of CD vs
#                control intestinal mucosa. Selects the soft-thresholding power
#                beta, builds signed co-expression modules, relates modules to
#                clinical traits (CD status, cohort) and characterises the
#                disease-associated "yellow" module through gene significance (GS)
#                vs module membership (MM).
#
# Inputs       : data/processed/GSE16879_expr_normalized.txt   (linear -> log2 here)
#                data/processed/GSE16879_phenotype.csv         (sample,group)
#                data/processed/GSE75214_expr_normalized.txt   (already log2)
#                data/processed/GSE75214_phenotype.csv
#                results/01_DEG/GSE16879_DEG_all.csv           (optional, for GS ranking)
#
# Outputs      : results/02_WGCNA/soft_threshold_table.csv
#                results/02_WGCNA/module_assignment.csv          (gene -> module colour)
#                results/02_WGCNA/module_trait_correlation.csv
#                results/02_WGCNA/module_eigengenes.csv
#                results/02_WGCNA/Yellow_Module_Genes_Info.csv   (GS, MM, hub flag)
#                results/02_WGCNA/hub_genes_yellow.csv
#                results/figures/WGCNA_SoftThreshold.{png,pdf}
#                results/figures/WGCNA_Gene_Dendrogram_Modules.{png,pdf}
#                results/figures/WGCNA_Module_Trait_Heatmap.{png,pdf}
#                results/figures/Yellow_Module_GS_vs_MM.{png,pdf}
#
# Figure/Table : Figure 2A-D ; Supplementary Table S3 (module membership)
#
# Key params   : MERGE_COHORTS = TRUE  -> discovery + validation merged on shared
#                                          genes with limma::removeBatchEffect(cohort)
#                N_TOP_GENES   = 5000  most variable genes (MAD)
#                networkType   = "signed" ; TOMType = "signed"
#                minModuleSize = 30 ; mergeCutHeight = 0.25 ; deepSplit = 2
#                beta selected as the smallest power with scale-free R^2 >= 0.85
#                (falls back to WGCNA::pickSoftThreshold powerEstimate)
#                MODULE_OF_INTEREST = "yellow" (CD-positively correlated module)
#                set.seed(2024)
#
# Runtime      : ~10-25 min (blockwiseModules on 5,000 genes, single block)
#                Memory: ~4 GB. TOM file is written next to the results folder.
#
# Packages     : WGCNA (>= 1.72), limma (>= 3.54), data.table (>= 1.14),
#                ggplot2 (>= 3.4), ggrepel (>= 0.9.3), flashClust/fastcluster
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

set.seed(2024)
options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()

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
RES_DIR  <- file.path(ROOT, "results", "02_WGCNA")
FIG_DIR  <- file.path(ROOT, "results", "figures")
DEG_DIR  <- file.path(ROOT, "results", "01_DEG")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

MERGE_COHORTS      <- TRUE
N_TOP_GENES        <- 5000
MIN_MODULE_SIZE    <- 30
MERGE_CUT_HEIGHT   <- 0.25
RSQ_CUT            <- 0.85
MODULE_OF_INTEREST <- "yellow"
HUB_GS             <- 0.2     # |GS| threshold for hub genes
HUB_MM             <- 0.8     # |MM| threshold for hub genes

## ---------------------------------------------------------------------------
## 1. Load and merge expression data
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
  mat
}

maybe_log2 <- function(mat, label = "") {
  mx <- max(mat, na.rm = TRUE)
  if (mx > 50) {
    message("  [", label, "] log2(x+1) applied (max was ", round(mx, 1), ")")
    mat[mat < 0] <- 0
    mat <- log2(mat + 1)
  }
  mat
}

read_phenotype <- function(path) {
  if (!file.exists(path)) stop("Missing phenotype file: ", path)
  ph <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  cn <- tolower(colnames(ph))
  scol <- colnames(ph)[which(cn %in% c("sample", "gsm", "geo_accession", "sample_id"))[1]]
  gcol <- colnames(ph)[which(cn %in% c("group", "class", "condition", "diagnosis", "disease"))[1]]
  data.frame(sample = as.character(ph[[scol]]),
             group  = ifelse(grepl("^(cd|crohn)", ph[[gcol]], ignore.case = TRUE), "CD",
                      ifelse(grepl("^(control|normal|non.?ibd|healthy)", ph[[gcol]],
                                   ignore.case = TRUE), "Control", NA)),
             stringsAsFactors = FALSE) |>
    subset(!is.na(group))
}

load_cohort <- function(expr_file, pheno_file, label) {
  e  <- maybe_log2(read_expr_matrix(file.path(DATA_DIR, expr_file)), label)
  ph <- read_phenotype(file.path(DATA_DIR, pheno_file))
  common <- intersect(colnames(e), ph$sample)
  ph <- ph[match(common, ph$sample), ]
  ph$cohort <- label
  list(expr = e[, common, drop = FALSE], pheno = ph)
}

message("=== Loading cohorts ===")
disc <- load_cohort("GSE16879_expr_normalized.txt", "GSE16879_phenotype.csv", "GSE16879")

if (MERGE_COHORTS && file.exists(file.path(DATA_DIR, "GSE75214_expr_normalized.txt"))) {
  vald   <- load_cohort("GSE75214_expr_normalized.txt", "GSE75214_phenotype.csv", "GSE75214")
  genes  <- intersect(rownames(disc$expr), rownames(vald$expr))
  message("  shared genes across cohorts: ", length(genes))
  # per-cohort quantile normalisation, then cross-cohort batch correction
  e1 <- limma::normalizeBetweenArrays(disc$expr[genes, , drop = FALSE], method = "quantile")
  e2 <- limma::normalizeBetweenArrays(vald$expr[genes, , drop = FALSE], method = "quantile")
  expr  <- cbind(e1, e2)
  pheno <- rbind(disc$pheno, vald$pheno)
  expr  <- limma::removeBatchEffect(
    expr, batch = factor(pheno$cohort),
    design = model.matrix(~ factor(pheno$group, levels = c("Control", "CD"))))
} else {
  expr  <- disc$expr
  pheno <- disc$pheno
}
stopifnot(identical(colnames(expr), pheno$sample))
message("  merged matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")

## ---------------------------------------------------------------------------
## 2. Gene / sample filtering
## ---------------------------------------------------------------------------
mads  <- apply(expr, 1, stats::mad, na.rm = TRUE)
sel   <- names(sort(mads, decreasing = TRUE))[seq_len(min(N_TOP_GENES, sum(mads > 0)))]
datExpr0 <- t(expr[sel, , drop = FALSE])          # WGCNA wants samples x genes

gsg <- WGCNA::goodSamplesGenes(datExpr0, verbose = 0)
if (!gsg$allOK) {
  message("  removing ", sum(!gsg$goodGenes), " genes / ", sum(!gsg$goodSamples), " samples (QC)")
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

# sample outlier detection by hierarchical clustering (documented, height-based)
sampleTree <- stats::hclust(stats::dist(datExpr0), method = "average")
cut_h      <- stats::quantile(sampleTree$height, 0.995) * 1.15
clust      <- stats::cutree(sampleTree, h = cut_h)
keep_samp  <- names(which(table(clust) == max(table(clust))))[1]
if (sum(clust == keep_samp) < nrow(datExpr0)) {
  message("  dropping ", sum(clust != keep_samp), " outlier sample(s) at h = ", round(cut_h, 1))
}
datExpr <- datExpr0[clust == keep_samp, , drop = FALSE]

ph <- pheno[match(rownames(datExpr), pheno$sample), ]
datTraits <- data.frame(
  CD        = as.numeric(ph$group == "CD"),
  Control   = as.numeric(ph$group == "Control"),
  row.names = ph$sample
)
if (length(unique(ph$cohort)) > 1) {
  datTraits$Validation_cohort <- as.numeric(ph$cohort == "GSE75214")
}
message("  WGCNA input: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

## ---------------------------------------------------------------------------
## 3. Soft-thresholding power (beta)
## ---------------------------------------------------------------------------
message("\n=== Picking soft-thresholding power ===")
powers <- c(1:10, seq(12, 30, by = 2))
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
                                networkType = "signed", verbose = 0)
sft_tab <- sft$fitIndices
sft_tab$scale_free_R2 <- -sign(sft_tab$slope) * sft_tab$SFT.R.sq
data.table::fwrite(sft_tab, file.path(RES_DIR, "soft_threshold_table.csv"))

ok    <- sft_tab$Power[sft_tab$scale_free_R2 >= RSQ_CUT]
BETA  <- if (length(ok)) min(ok) else if (!is.na(sft$powerEstimate)) sft$powerEstimate else 12
message("  selected beta = ", BETA,
        " (scale-free R^2 = ", round(sft_tab$scale_free_R2[sft_tab$Power == BETA], 3), ")")

png(file.path(FIG_DIR, "WGCNA_SoftThreshold.png"), width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(sft_tab$Power, sft_tab$scale_free_R2, type = "n",
     xlab = "Soft-threshold power (beta)", ylab = "Scale-free topology model fit, signed R^2",
     main = "Scale independence")
text(sft_tab$Power, sft_tab$scale_free_R2, labels = sft_tab$Power, col = "#E94F37", cex = 0.9)
abline(h = RSQ_CUT, col = "#2E86AB", lty = 2)
plot(sft_tab$Power, sft_tab$mean.k., type = "n",
     xlab = "Soft-threshold power (beta)", ylab = "Mean connectivity", main = "Mean connectivity")
text(sft_tab$Power, sft_tab$mean.k., labels = sft_tab$Power, col = "#E94F37", cex = 0.9)
dev.off()
pdf(file.path(FIG_DIR, "WGCNA_SoftThreshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(sft_tab$Power, sft_tab$scale_free_R2, type = "n", xlab = "Soft-threshold power (beta)",
     ylab = "Scale-free topology model fit, signed R^2", main = "Scale independence")
text(sft_tab$Power, sft_tab$scale_free_R2, labels = sft_tab$Power, col = "#E94F37", cex = 0.9)
abline(h = RSQ_CUT, col = "#2E86AB", lty = 2)
plot(sft_tab$Power, sft_tab$mean.k., type = "n", xlab = "Soft-threshold power (beta)",
     ylab = "Mean connectivity", main = "Mean connectivity")
text(sft_tab$Power, sft_tab$mean.k., labels = sft_tab$Power, col = "#E94F37", cex = 0.9)
dev.off()

## ---------------------------------------------------------------------------
## 4. Network construction and module detection
## ---------------------------------------------------------------------------
message("\n=== Building network (blockwiseModules) ===")
net <- WGCNA::blockwiseModules(
  datExpr,
  power              = BETA,
  networkType        = "signed",
  TOMType            = "signed",
  maxBlockSize       = ncol(datExpr),
  minModuleSize      = MIN_MODULE_SIZE,
  deepSplit          = 2,
  reassignThreshold  = 0,
  mergeCutHeight     = MERGE_CUT_HEIGHT,
  numericLabels      = TRUE,
  pamRespectsDendro  = FALSE,
  saveTOMs           = TRUE,
  saveTOMFileBase    = file.path(RES_DIR, "CD_vs_Control_TOM"),
  randomSeed         = 2024,
  verbose            = 3
)

moduleColors <- WGCNA::labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
message("  modules detected: ", length(unique(moduleColors)))
print(sort(table(moduleColors), decreasing = TRUE))

data.table::fwrite(
  data.frame(Gene = names(moduleColors), Module = as.character(moduleColors),
             ModuleLabel = net$colors[names(moduleColors)]),
  file.path(RES_DIR, "module_assignment.csv"))

# Dendrogram + module colours
png(file.path(FIG_DIR, "WGCNA_Gene_Dendrogram_Modules.png"),
    width = 2200, height = 1300, res = 200)
WGCNA::plotDendroAndColors(net$dendrograms[[1]],
                           moduleColors[net$blockGenes[[1]]],
                           "Module colours", dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05,
                           main = "Gene clustering dendrogram and module colours")
dev.off()
pdf(file.path(FIG_DIR, "WGCNA_Gene_Dendrogram_Modules.pdf"), width = 11, height = 6.5)
WGCNA::plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                           "Module colours", dendroLabels = FALSE, hang = 0.03,
                           addGuide = TRUE, guideHang = 0.05,
                           main = "Gene clustering dendrogram and module colours")
dev.off()

## ---------------------------------------------------------------------------
## 5. Module-trait relationships
## ---------------------------------------------------------------------------
message("\n=== Module-trait relationships ===")
MEs0 <- WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs  <- WGCNA::orderMEs(MEs0)
data.table::fwrite(data.frame(sample = rownames(datExpr), MEs),
                   file.path(RES_DIR, "module_eigengenes.csv"))

moduleTraitCor  <- stats::cor(MEs, datTraits, use = "pairwise.complete.obs")
moduleTraitP    <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(datExpr))

mt <- data.frame(Module = rownames(moduleTraitCor), moduleTraitCor, check.names = FALSE)
colnames(mt)[-1] <- paste0("cor_", colnames(moduleTraitCor))
mtp <- as.data.frame(moduleTraitP); colnames(mtp) <- paste0("p_", colnames(moduleTraitP))
data.table::fwrite(cbind(mt, mtp), file.path(RES_DIR, "module_trait_correlation.csv"))

textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitP, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

draw_mt_heatmap <- function() {
  par(mar = c(6, 9, 3, 1.5))
  WGCNA::labeledHeatmap(
    Matrix = moduleTraitCor, xLabels = colnames(datTraits), yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor), colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50), textMatrix = textMatrix,
    setStdMargins = FALSE, cex.text = 0.65, zlim = c(-1, 1),
    main = "Module-trait relationships")
}
png(file.path(FIG_DIR, "WGCNA_Module_Trait_Heatmap.png"), width = 1700, height = 2000, res = 200)
draw_mt_heatmap(); dev.off()
pdf(file.path(FIG_DIR, "WGCNA_Module_Trait_Heatmap.pdf"), width = 8, height = 10)
draw_mt_heatmap(); dev.off()

cd_cor <- moduleTraitCor[, "CD"]
message("  top CD-correlated modules:")
print(round(head(sort(cd_cor, decreasing = TRUE), 5), 3))

## ---------------------------------------------------------------------------
## 6. Gene significance (GS) vs module membership (MM) -- yellow module
## ---------------------------------------------------------------------------
message("\n=== GS vs MM for the '", MODULE_OF_INTEREST, "' module ===")
modNames <- substring(colnames(MEs), 3)

geneModuleMembership <- as.data.frame(stats::cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
colnames(geneModuleMembership) <- paste0("MM", modNames)
colnames(MMPvalue)             <- paste0("p.MM", modNames)

trait_CD <- as.data.frame(datTraits$CD); names(trait_CD) <- "CD"
geneTraitSignificance <- as.data.frame(stats::cor(datExpr, trait_CD, use = "p"))
GSPvalue <- as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))
colnames(geneTraitSignificance) <- "GS.CD"
colnames(GSPvalue)              <- "p.GS.CD"

target_module <- if (MODULE_OF_INTEREST %in% modNames) {
  MODULE_OF_INTEREST
} else {
  warning("Module '", MODULE_OF_INTEREST, "' not found; using the most CD-correlated module.")
  substring(names(which.max(cd_cor)), 3)
}
column   <- match(target_module, modNames)
inModule <- moduleColors == target_module
message("  genes in '", target_module, "' module: ", sum(inModule))

mod_info <- data.frame(
  Gene   = colnames(datExpr)[inModule],
  Module = target_module,
  MM     = geneModuleMembership[inModule, column],
  MM_p   = MMPvalue[inModule, column],
  GS     = geneTraitSignificance$GS.CD[inModule],
  GS_p   = GSPvalue$p.GS.CD[inModule],
  stringsAsFactors = FALSE
)
mod_info$absMM <- abs(mod_info$MM); mod_info$absGS <- abs(mod_info$GS)
mod_info$Hub   <- mod_info$absMM > HUB_MM & mod_info$absGS > HUB_GS
mod_info <- mod_info[order(-mod_info$absMM, -mod_info$absGS), ]

# add limma statistics of the discovery cohort when available
deg_file <- file.path(DEG_DIR, "GSE16879_DEG_all.csv")
if (file.exists(deg_file)) {
  deg <- data.table::fread(deg_file, data.table = FALSE)
  mod_info$logFC_GSE16879 <- deg$logFC[match(mod_info$Gene, deg$Gene)]
  mod_info$adjP_GSE16879  <- deg$adj.P.Val[match(mod_info$Gene, deg$Gene)]
}
data.table::fwrite(mod_info, file.path(RES_DIR,
  paste0(tools::toTitleCase(target_module), "_Module_Genes_Info.csv")))
data.table::fwrite(mod_info[mod_info$Hub, ],
                   file.path(RES_DIR, paste0("hub_genes_", target_module, ".csv")))
message("  hub genes (|MM| > ", HUB_MM, " & |GS| > ", HUB_GS, "): ", sum(mod_info$Hub))

cor_gsmm <- stats::cor.test(mod_info$MM, mod_info$GS)
label_genes <- unique(c(
  head(mod_info$Gene[mod_info$Hub], 12),
  intersect(c("CFI", "PAQR5", "KCNE3", "C2", "CFB"), mod_info$Gene)))

p_gsmm <- ggplot(mod_info, aes(x = MM, y = GS)) +
  geom_point(colour = target_module, alpha = 0.65, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey25", linewidth = 0.6, formula = y ~ x) +
  ggrepel::geom_text_repel(data = subset(mod_info, Gene %in% label_genes),
                           aes(label = Gene), size = 3, max.overlaps = 40,
                           box.padding = 0.35, segment.size = 0.25) +
  labs(title = paste0("Module membership vs gene significance (", target_module, " module)"),
       subtitle = sprintf("Pearson r = %.2f, P = %.2g, n = %d genes",
                          cor_gsmm$estimate, cor_gsmm$p.value, nrow(mod_info)),
       x = paste0("Module membership (MM) in the ", target_module, " module"),
       y = "Gene significance for CD") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Yellow_Module_GS_vs_MM.png"), p_gsmm, width = 7.5, height = 6.5, dpi = 300)
ggsave(file.path(FIG_DIR, "Yellow_Module_GS_vs_MM.pdf"), p_gsmm, width = 7.5, height = 6.5)

message("\nDone. Tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)

## ---------------------------------------------------------------------------
## 7. Session information
## ---------------------------------------------------------------------------
sessionInfo()
