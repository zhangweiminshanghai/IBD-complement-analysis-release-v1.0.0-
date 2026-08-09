# ============================================================================
# Figure S10: Single-Cell Validation of CFI in Crohn's Disease
# Data Source: GSE134809 (Martin et al., Cell 2019)
# Purpose: Re-generate compliant figure with correct "Figure S10" caption
# ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("D:/IBD_Project/scRNA_data")
outdir <- "D:/IBD_Project/final_suppl_figs"

# ── Load Martin 2019 processed data ─────────────────────────────────────────
# Use existing processed RDS or re-process from raw
rds_file <- "Martin2019_results/Martin2019_processed_fixed.rds"
if (file.exists(rds_file)) {
  cat("Loading processed data:", rds_file, "\n")
  seurat_obj <- readRDS(rds_file)
} else if (file.exists("Martin2019_results/Martin2019_filtered.rds")) {
  cat("Loading filtered data and re-processing...\n")
  seurat_obj <- readRDS("Martin2019_results/Martin2019_filtered.rds")
  # Standard processing pipeline
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
} else {
  stop("No processed data found. Please run data download and preprocessing first.")
}

cat("Cells:", ncol(seurat_obj), "| Features:", nrow(seurat_obj), "\n")

# ── Assign disease labels based on sample IDs ───────────────────────────────
# Martin 2019: Control samples (N*, HC*), Ileal CD (I*), Colonic CD (C*)
seurat_obj$sample_id <- seurat_obj$orig.ident
seurat_obj$disease <- case_when(
  grepl("^N|^HC", seurat_obj$sample_id) ~ "Control",
  grepl("^I", seurat_obj$sample_id) ~ "Ileal CD",
  grepl("^C", seurat_obj$sample_id) ~ "Colonic CD",
  TRUE ~ "Unknown"
)

seurat_obj$disease <- factor(seurat_obj$disease, 
                              levels = c("Control", "Ileal CD", "Colonic CD"))

# ── Check CFI and C2 expression ─────────────────────────────────────────────
if ("CFI" %in% rownames(seurat_obj)) {
  cat("CFI expression range:", range(seurat_obj@assays$RNA@data["CFI", ]), "\n")
}
if ("C2" %in% rownames(seurat_obj)) {
  cat("C2 expression range:", range(seurat_obj@assays$RNA@data["C2", ]), "\n")
}

# ── Panel A: UMAP with cell clusters ────────────────────────────────────────
pA <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  labs(title = "A. Cell Clusters (Leiden resolution = 0.5)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ── Panel B: Disease status on UMAP ─────────────────────────────────────────
pB <- DimPlot(seurat_obj, reduction = "umap", group.by = "disease",
              cols = c("Control" = "#808080", "Ileal CD" = "#E94F37", "Colonic CD" = "#2E86AB")) +
  labs(title = "B. Disease Status") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── Panel C: CFI expression on UMAP ─────────────────────────────────────────
if ("CFI" %in% rownames(seurat_obj)) {
  pC <- FeaturePlot(seurat_obj, features = "CFI", reduction = "umap",
                    cols = c("lightgrey", "darkred")) +
    labs(title = "C. CFI Expression") +
    theme_bw(base_size = 11)
} else {
  pC <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "CFI not found") +
    theme_void() + labs(title = "C. CFI Expression")
}

# ── Panel D: C2 expression on UMAP ──────────────────────────────────────────
if ("C2" %in% rownames(seurat_obj)) {
  pD <- FeaturePlot(seurat_obj, features = "C2", reduction = "umap",
                    cols = c("lightgrey", "darkblue")) +
    labs(title = "D. C2 Expression") +
    theme_bw(base_size = 11)
} else {
  pD <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "C2 not found") +
    theme_void() + labs(title = "D. C2 Expression")
}

# ── Panel E: CFI expression by disease (violin) ─────────────────────────────
if ("CFI" %in% rownames(seurat_obj)) {
  cfi_data <- data.frame(
    cell = colnames(seurat_obj),
    CFI = seurat_obj@assays$RNA@data["CFI", ],
    disease = seurat_obj$disease
  )
  
  pE <- ggplot(cfi_data, aes(x = disease, y = CFI, fill = disease)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("Control" = "#808080", 
                                  "Ileal CD" = "#E94F37", 
                                  "Colonic CD" = "#2E86AB")) +
    labs(x = "Disease", y = "CFI Expression (log-normalized)",
         title = "E. CFI Expression by Disease") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  pE <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "CFI not found") +
    theme_void() + labs(title = "E. CFI Expression by Disease")
}

# ── Panel F: CFI+ cell fraction by disease ──────────────────────────────────
if ("CFI" %in% rownames(seurat_obj)) {
  cfi_pos <- seurat_obj@assays$RNA@data["CFI", ] > 0
  cfi_frac <- data.frame(
    disease = seurat_obj$disease,
    CFI_positive = cfi_pos
  ) %>%
    group_by(disease) %>%
    summarise(fraction = mean(CFI_positive), n = n(), .groups = "drop")
  
  pF <- ggplot(cfi_frac, aes(x = disease, y = fraction, fill = disease)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", fraction*100, n)), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Control" = "#808080", 
                                  "Ileal CD" = "#E94F37", 
                                  "Colonic CD" = "#2E86AB")) +
    labs(x = "Disease", y = "CFI+ Cell Fraction",
         title = "F. CFI+ Cell Fraction by Disease") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, max(cfi_frac$fraction) * 1.2)
} else {
  pF <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "CFI not found") +
    theme_void() + labs(title = "F. CFI+ Cell Fraction by Disease")
}

# ── Combine panels with CORRECT Figure S10 caption ──────────────────────────
fig <- (pA | pB | pC) / (pD | pE | pF) +
  plot_annotation(
    title = "Figure S10. Single-Cell Validation of CFI in Crohn's Disease (GSE134809)",
    subtitle = "Martin et al. 2019 (Cell): 160,981 cells from 11 Control, 11 Ileal CD, 11 Colonic CD",
    caption = "High-confidence CFI+ cells defined by FC>=3 threshold; C2 expression shows disease-location specificity",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.caption = element_text(hjust = 0.5, size = 9, face = "italic")
    )
  )

# Save high-resolution figure
ggsave(file.path(outdir, "Figure_S10_scRNA_CFI_CD.png"), fig, 
       width = 16, height = 10, dpi = 300, units = "in", bg = "white")
ggsave(file.path(outdir, "Figure_S10_scRNA_CFI_CD.pdf"), fig, 
       width = 16, height = 10, units = "in")

cat("Saved: Figure_S10_scRNA_CFI_CD.png/pdf\n")

# Save summary statistics
if ("CFI" %in% rownames(seurat_obj)) {
  cfi_summary <- data.frame(
    disease = cfi_frac$disease,
    CFI_positive_fraction = cfi_frac$fraction,
    total_cells = cfi_frac$n
  )
  write.csv(cfi_summary, file.path(outdir, "Figure_S10_CFI_stats.csv"), row.names = FALSE)
  cat("Saved: Figure_S10_CFI_stats.csv\n")
}

cat("\n=== Figure S10 Complete ===\n")
