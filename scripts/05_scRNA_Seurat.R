###############################################################################
# 05_scRNA_Seurat.R
#
# Purpose      : Single-cell validation of CFI expression in Crohn's disease
#                (GSE134809, Martin et al. Cell 2019; 160,981 cells, 22 samples).
#                Full Seurat v5 pipeline: load 10x matrices, QC, normalise,
#                integrate, cluster (23 clusters), UMAP, then define CFI+ cells
#                with two thresholds and test enrichment in ileal CD.
#
# Inputs       : data/processed/GSE134809/                      (one sub-directory per
#                    sample containing barcodes.tsv(.gz), features/genes.tsv(.gz),
#                    matrix.mtx(.gz)), OR
#                data/processed/GSE134809_seurat.rds            (pre-built Seurat object)
#                data/processed/GSE134809_sample_metadata.csv   (optional; columns:
#                    sample,patient,tissue,disease  with disease in
#                    {Control, Ileal CD, Colonic CD})
#
# Outputs      : results/05_scRNA/GSE134809_seurat_processed.rds
#                results/05_scRNA/cluster_composition.csv
#                results/05_scRNA/CFIpos_cells_by_group.csv
#                results/05_scRNA/CFIpos_cells_by_cluster.csv
#                results/05_scRNA/CFIpos_enrichment_tests.csv
#                results/05_scRNA/cell_metadata.csv.gz
#                results/figures/scRNA_UMAP_clusters.{png,pdf}
#                results/figures/scRNA_UMAP_disease.{png,pdf}
#                results/figures/scRNA_UMAP_CFI.{png,pdf}
#                results/figures/scRNA_CFIpos_proportions.{png,pdf}
#
# Figure/Table : Figure 5A-D ; Supplementary Table S6 (CFI+ proportions)
#
# Key params   : QC        : 200 <= nFeature_RNA <= 6000, nCount_RNA >= 500,
#                            percent.mt < 20
#                Normalise : LogNormalize, scale.factor = 1e4
#                HVG       : vst, nfeatures = 2000
#                PCA/UMAP  : 30 PCs ; integration = Harmony (per sample) when available
#                Clustering: Louvain; resolution auto-tuned to TARGET_CLUSTERS = 23
#                CFI+ std  : raw UMI count of CFI > 0
#                CFI+ high : raw UMI count >= 90th percentile of CFI+ cells (top decile)
#                Enrichment: two-sided Fisher exact test, ileal CD vs all other cells
#                            (reported: p = 2.4e-22; 0.89% of ileal CD cells CFI+;
#                             CFI+ cells absent from colonic CD)
#                set.seed(2024)
#
# Runtime      : ~60-120 min and >= 64 GB RAM for the full 160k-cell object
#                (~10 min if a pre-built RDS is supplied).
#
# Packages     : Seurat (>= 5.0.1), SeuratObject (>= 5.0), Matrix (>= 1.6),
#                dplyr (>= 1.1), ggplot2 (>= 3.4), patchwork (>= 1.1),
#                data.table (>= 1.14); optional harmony (>= 1.2), presto
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

set.seed(2024)
options(future.globals.maxSize = 8 * 1024^3)

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
RES_DIR  <- file.path(ROOT, "results", "05_scRNA")
FIG_DIR  <- file.path(ROOT, "results", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SC_DIR          <- file.path(DATA_DIR, "GSE134809")
RDS_IN          <- file.path(DATA_DIR, "GSE134809_seurat.rds")
META_FILE       <- file.path(DATA_DIR, "GSE134809_sample_metadata.csv")
RDS_OUT         <- file.path(RES_DIR, "GSE134809_seurat_processed.rds")

MIN_FEATURES    <- 200
MAX_FEATURES    <- 6000
MIN_COUNTS      <- 500
MAX_MT          <- 20
N_HVG           <- 2000
N_PCS           <- 30
TARGET_CLUSTERS <- 23
RES_GRID        <- seq(0.2, 2.0, by = 0.1)
GENE            <- "CFI"
TOP_DECILE      <- 0.90

## ---------------------------------------------------------------------------
## 1. Load data
## ---------------------------------------------------------------------------
load_object <- function() {
  if (file.exists(RDS_IN)) {
    message("Loading pre-built Seurat object: ", RDS_IN)
    return(readRDS(RDS_IN))
  }
  if (!dir.exists(SC_DIR)) {
    stop("Neither ", RDS_IN, " nor ", SC_DIR, " exists.\n",
         "  Download GSE134809 (Martin et al. 2019) 10x matrices into\n",
         "  data/processed/GSE134809/<sample>/{barcodes,features,matrix}.")
  }
  samples <- list.dirs(SC_DIR, recursive = FALSE, full.names = TRUE)
  message("Reading ", length(samples), " 10x sample directories ...")
  obj_list <- lapply(samples, function(sdir) {
    sid <- basename(sdir)
    cnt <- Seurat::Read10X(sdir)
    if (is.list(cnt)) cnt <- cnt[[1]]
    o <- Seurat::CreateSeuratObject(counts = cnt, project = sid,
                                    min.cells = 3, min.features = 100)
    o$sample_id <- sid
    o
  })
  if (length(obj_list) == 1) return(obj_list[[1]])
  merged <- merge(obj_list[[1]], y = obj_list[-1],
                  add.cell.ids = basename(samples), project = "GSE134809")
  merged
}

seu <- load_object()
message("Raw object: ", ncol(seu), " cells x ", nrow(seu), " features")

if (is.null(seu$sample_id)) seu$sample_id <- as.character(Seurat::Idents(seu))

## ---------------------------------------------------------------------------
## 2. Sample-level metadata (tissue / disease group)
## ---------------------------------------------------------------------------
assign_disease <- function(sample_ids) {
  # Martin 2019 naming: N*/HC* = control, I* = ileal CD, C* = colonic CD
  dplyr::case_when(
    grepl("colon", sample_ids, ignore.case = TRUE)  ~ "Colonic CD",
    grepl("ile",   sample_ids, ignore.case = TRUE)  ~ "Ileal CD",
    grepl("^(N|HC|Ctrl|Control)", sample_ids)       ~ "Control",
    grepl("^I", sample_ids)                         ~ "Ileal CD",
    grepl("^C", sample_ids)                         ~ "Colonic CD",
    TRUE                                            ~ "Unknown")
}

if (file.exists(META_FILE)) {
  meta <- data.table::fread(META_FILE, data.table = FALSE)
  idx  <- match(seu$sample_id, meta$sample)
  seu$patient <- if ("patient" %in% colnames(meta)) meta$patient[idx] else seu$sample_id
  seu$tissue  <- if ("tissue"  %in% colnames(meta)) meta$tissue[idx]  else NA_character_
  seu$disease <- if ("disease" %in% colnames(meta)) meta$disease[idx] else assign_disease(seu$sample_id)
} else {
  message("No sample metadata CSV found -> deriving disease labels from sample IDs.")
  seu$patient <- seu$sample_id
  seu$disease <- assign_disease(seu$sample_id)
}
seu$disease <- factor(seu$disease,
                      levels = intersect(c("Control", "Ileal CD", "Colonic CD", "Unknown"),
                                         unique(as.character(seu$disease))))
print(table(seu$disease, useNA = "ifany"))

## ---------------------------------------------------------------------------
## 3. Quality control
## ---------------------------------------------------------------------------
seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
qc_before <- ncol(seu)
seu <- subset(seu, subset = nFeature_RNA >= MIN_FEATURES & nFeature_RNA <= MAX_FEATURES &
                            nCount_RNA  >= MIN_COUNTS   & percent.mt   <  MAX_MT)
message("QC: ", qc_before, " -> ", ncol(seu), " cells retained")

p_qc <- Seurat::VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        group.by = "sample_id", pt.size = 0, ncol = 3) &
  theme(axis.text.x = element_text(size = 6, angle = 90))
ggsave(file.path(FIG_DIR, "scRNA_QC_violin.png"), p_qc, width = 16, height = 5, dpi = 200)

## ---------------------------------------------------------------------------
## 4. Normalisation, HVG, scaling, PCA (Seurat v5 layer API)
## ---------------------------------------------------------------------------
seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = N_HVG)
seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu),
                         vars.to.regress = NULL)
seu <- Seurat::RunPCA(seu, features = Seurat::VariableFeatures(seu),
                      npcs = N_PCS, seed.use = 2024, verbose = FALSE)

reduction_use <- "pca"
if (length(unique(seu$sample_id)) > 1 && requireNamespace("harmony", quietly = TRUE)) {
  message("Running Harmony batch integration over sample_id ...")
  seu <- harmony::RunHarmony(seu, group.by.vars = "sample_id",
                             reduction.use = "pca", dims.use = 1:N_PCS,
                             plot_convergence = FALSE)
  reduction_use <- "harmony"
} else if (inherits(seu[["RNA"]], "Assay5") &&
           length(SeuratObject::Layers(seu, search = "counts")) > 1) {
  message("Running Seurat v5 CCA integration ...")
  seu <- Seurat::IntegrateLayers(seu, method = Seurat::CCAIntegration,
                                 orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)
  reduction_use <- "integrated.cca"
}

## ---------------------------------------------------------------------------
## 5. Clustering tuned to 23 clusters + UMAP
## ---------------------------------------------------------------------------
seu <- Seurat::FindNeighbors(seu, reduction = reduction_use, dims = 1:N_PCS, verbose = FALSE)

chosen_res <- NA_real_
for (r in RES_GRID) {
  seu <- Seurat::FindClusters(seu, resolution = r, random.seed = 2024, verbose = FALSE)
  k <- length(unique(Seurat::Idents(seu)))
  message(sprintf("  resolution %.2f -> %d clusters", r, k))
  chosen_res <- r
  if (k >= TARGET_CLUSTERS) break
}
message("Selected resolution = ", chosen_res, " with ",
        length(unique(Seurat::Idents(seu))), " clusters (target ", TARGET_CLUSTERS, ")")
seu$seurat_clusters <- Seurat::Idents(seu)

seu <- Seurat::RunUMAP(seu, reduction = reduction_use, dims = 1:N_PCS,
                       seed.use = 2024, verbose = FALSE)

if (inherits(seu[["RNA"]], "Assay5")) seu <- SeuratObject::JoinLayers(seu)

## ---------------------------------------------------------------------------
## 6. CFI+ cell definitions
## ---------------------------------------------------------------------------
if (!GENE %in% rownames(seu)) stop(GENE, " is not present in the expression matrix.")
counts_cfi <- as.numeric(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")[GENE, ])
norm_cfi   <- as.numeric(SeuratObject::LayerData(seu, assay = "RNA", layer = "data")[GENE, ])

seu$CFI_counts <- counts_cfi
seu$CFI_norm   <- norm_cfi

# (a) standard definition: any raw UMI
seu$CFI_pos <- counts_cfi > 0
# (b) high-confidence definition: top decile among CFI-expressing cells
pos_vals <- counts_cfi[counts_cfi > 0]
cut_high <- if (length(pos_vals)) as.numeric(stats::quantile(pos_vals, TOP_DECILE)) else Inf
seu$CFI_pos_high <- counts_cfi >= cut_high & counts_cfi > 0

message(sprintf("CFI+ (UMI > 0): %d cells (%.2f%%)",
                sum(seu$CFI_pos), 100 * mean(seu$CFI_pos)))
message(sprintf("CFI+ high-confidence (>= %g UMI, top decile): %d cells (%.2f%%)",
                cut_high, sum(seu$CFI_pos_high), 100 * mean(seu$CFI_pos_high)))

## ---------------------------------------------------------------------------
## 7. Proportion tables and enrichment tests
## ---------------------------------------------------------------------------
md <- seu@meta.data
by_group <- md %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(n_cells       = dplyr::n(),
                   n_CFIpos      = sum(CFI_pos),
                   pct_CFIpos    = 100 * mean(CFI_pos),
                   n_CFIpos_high = sum(CFI_pos_high),
                   pct_CFIpos_high = 100 * mean(CFI_pos_high),
                   .groups = "drop")
data.table::fwrite(by_group, file.path(RES_DIR, "CFIpos_cells_by_group.csv"))
print(as.data.frame(by_group))

by_cluster <- md %>%
  dplyr::group_by(seurat_clusters, disease) %>%
  dplyr::summarise(n_cells = dplyr::n(), n_CFIpos = sum(CFI_pos),
                   pct_CFIpos = 100 * mean(CFI_pos),
                   mean_CFI_norm = mean(CFI_norm), .groups = "drop")
data.table::fwrite(by_cluster, file.path(RES_DIR, "CFIpos_cells_by_cluster.csv"))

comp <- md %>% dplyr::count(seurat_clusters, disease, name = "n_cells") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(pct_of_cluster = 100 * n_cells / sum(n_cells)) %>% dplyr::ungroup()
data.table::fwrite(comp, file.path(RES_DIR, "cluster_composition.csv"))

# Fisher exact test: CFI+ enrichment of each disease group vs all other cells
test_rows <- lapply(levels(droplevels(factor(md$disease))), function(grp) {
  tab <- table(factor(md$disease == grp, levels = c(FALSE, TRUE)),
               factor(md$CFI_pos,        levels = c(FALSE, TRUE)))
  if (any(dim(tab) < 2)) return(NULL)
  ft <- stats::fisher.test(tab)
  ct <- suppressWarnings(stats::chisq.test(tab))
  data.frame(Group = grp,
             n_cells = sum(md$disease == grp),
             n_CFIpos = sum(md$disease == grp & md$CFI_pos),
             pct_CFIpos = 100 * mean(md$CFI_pos[md$disease == grp]),
             odds_ratio = unname(ft$estimate),
             fisher_p = ft$p.value, chisq_p = ct$p.value,
             stringsAsFactors = FALSE)
})
tests <- do.call(rbind, Filter(Negate(is.null), test_rows))
tests$fisher_FDR <- stats::p.adjust(tests$fisher_p, method = "BH")
data.table::fwrite(tests, file.path(RES_DIR, "CFIpos_enrichment_tests.csv"))
print(tests)

if ("Ileal CD" %in% tests$Group) {
  i <- which(tests$Group == "Ileal CD")
  message(sprintf("Ileal CD: %.2f%% of cells are CFI+ (Fisher P = %.3g, OR = %.2f)",
                  tests$pct_CFIpos[i], tests$fisher_p[i], tests$odds_ratio[i]))
  message("Reference values from the manuscript: 0.89% of ileal CD cells, p = 2.4e-22.")
}
if ("Colonic CD" %in% tests$Group) {
  i <- which(tests$Group == "Colonic CD")
  message(sprintf("Colonic CD: %d CFI+ cells (%.3f%%) - expected to be absent.",
                  tests$n_CFIpos[i], tests$pct_CFIpos[i]))
}

data.table::fwrite(data.frame(cell = rownames(md), md),
                   file.path(RES_DIR, "cell_metadata.csv.gz"))

## ---------------------------------------------------------------------------
## 8. Figures
## ---------------------------------------------------------------------------
p_clusters <- Seurat::DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE,
                              raster = FALSE) +
  labs(title = sprintf("GSE134809: %s cells, %d clusters (res = %.1f)",
                       format(ncol(seu), big.mark = ","),
                       length(unique(seu$seurat_clusters)), chosen_res)) +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_disease <- Seurat::DimPlot(seu, reduction = "umap", group.by = "disease", raster = FALSE,
                             cols = c("Control" = "#8A8A8A", "Ileal CD" = "#E94F37",
                                      "Colonic CD" = "#2E86AB", "Unknown" = "#CCCCCC")) +
  labs(title = "Disease group") + theme_bw(base_size = 11)

p_cfi <- Seurat::FeaturePlot(seu, features = GENE, reduction = "umap", raster = FALSE,
                             cols = c("grey88", "#B3122E"), order = TRUE) +
  labs(title = paste0(GENE, " expression (log-normalised)")) + theme_bw(base_size = 11)

for (nm in c("clusters", "disease", "CFI")) {
  p <- switch(nm, clusters = p_clusters, disease = p_disease, CFI = p_cfi)
  ggsave(file.path(FIG_DIR, paste0("scRNA_UMAP_", nm, ".png")), p,
         width = 8, height = 7, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0("scRNA_UMAP_", nm, ".pdf")), p, width = 8, height = 7)
}

p_prop <- ggplot(by_group, aes(x = disease, y = pct_CFIpos, fill = disease)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f%%\n(n = %d)", pct_CFIpos, n_CFIpos)),
            vjust = -0.25, size = 3.2) +
  scale_fill_manual(values = c("Control" = "#8A8A8A", "Ileal CD" = "#E94F37",
                               "Colonic CD" = "#2E86AB", "Unknown" = "#CCCCCC"),
                    guide = "none") +
  expand_limits(y = max(by_group$pct_CFIpos) * 1.25) +
  labs(title = "CFI+ cells by disease group",
       subtitle = "CFI+ defined as raw UMI count > 0 (GSE134809)",
       x = NULL, y = "% CFI+ cells") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "scRNA_CFIpos_proportions.png"), p_prop,
       width = 6.5, height = 5.5, dpi = 300)
ggsave(file.path(FIG_DIR, "scRNA_CFIpos_proportions.pdf"), p_prop, width = 6.5, height = 5.5)

combined <- (p_clusters | p_disease) / (p_cfi | p_prop) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave(file.path(FIG_DIR, "scRNA_Figure5_composite.png"), combined,
       width = 15, height = 13, dpi = 300)
ggsave(file.path(FIG_DIR, "scRNA_Figure5_composite.pdf"), combined, width = 15, height = 13)

saveRDS(seu, RDS_OUT)
message("\nDone. Object -> ", RDS_OUT, " ; tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)

## ---------------------------------------------------------------------------
## 9. Session information
## ---------------------------------------------------------------------------
sessionInfo()
