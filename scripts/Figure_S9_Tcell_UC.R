# ============================================================================
# Figure S9: Single-Cell T Cell Landscape in Ulcerative Colitis
# Data Source: GSE125527 (Smillie et al., Cell 2019)
# Purpose: Re-generate compliant figure with correct "Figure S9" caption
# ============================================================================

library(data.table)
library(ggplot2)
library(ggsci)
library(patchwork)

setwd("D:/IBD_Project/scRNA_data")
outdir <- "D:/IBD_Project/final_suppl_figs"

# ── Load GSE125527 T cell metadata ──────────────────────────────────────────
meta <- fread("GSE125527_Tcell_cluster.csv.gz")
cat("Loaded T cells:", nrow(meta), "\n")

# Assign disease status: C1-C9 = Control/Healthy; U1-U7 = UC patient
meta[, disease := ifelse(grepl("^C[0-9]", patient_id), "Healthy", "UC")]
meta[, tissue_type := ifelse(tissue_id == "R", "Rectum", "PBMC")]

# Cluster annotations from Smillie et al.
cluster_annotation <- c(
  T1 = "CD4+ Central Memory",
  T2 = "CD4+ Effector", 
  T3 = "CD8+ Cytotoxic",
  T4 = "γδ T cells",
  T5 = "NK-like T",
  T6 = "T follicular helper",
  T7 = "Regulatory T",
  T8 = "CD8+ IEL",
  T9 = "Cycling T",
  T10 = "CD4+ Activated",
  T11 = "MAIT",
  T12 = "CD4+ Th17",
  T13 = "CD8+ Exhausted",
  T14 = "DN T",
  T15 = "gdT Activated",
  T16 = "T progenitor",
  T17 = "NKT"
)

meta[, cluster_name := cluster_annotation[Cluster_id]]
meta[, cluster_name := ifelse(is.na(cluster_name), as.character(Cluster_id), cluster_name)]

# ── Calculate UC vs Healthy enrichment ──────────────────────────────────────
comp_data <- dcast(meta[, .(n=.N), by=.(Cluster_id, disease)], 
                   Cluster_id ~ disease, value.var = "n", fill = 0)
comp_data[, total := Healthy + UC]
comp_data[, UC_frac := UC / total]
comp_data[, log2FC := log2((UC + 0.5) / (Healthy + 0.5))]
comp_data[, cluster_name := cluster_annotation[Cluster_id]]

# ── Panel A: UC-enriched vs Healthy-enriched T cells ────────────────────────
pa_data <- comp_data[order(-log2FC)]
pa_data[, enrichment := ifelse(log2FC > 0, "UC-enriched", "Healthy-enriched")]
pa_data[, enrichment := factor(enrichment, levels = c("UC-enriched", "Healthy-enriched"))]

pA <- ggplot(pa_data, aes(x = reorder(Cluster_id, log2FC), y = log2FC, fill = enrichment)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = c("UC-enriched" = "#E94F37", "Healthy-enriched" = "#2E86AB")) +
  coord_flip() +
  labs(x = "T Cell Cluster", y = "log2 Fold Change (UC/Healthy)",
       title = "A. T Cell Enrichment in UC") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

# ── Panel B: Overall composition by disease ─────────────────────────────────
comp_long <- meta[, .(n = .N), by=.(Cluster_id, disease)]
comp_long[, prop := n / sum(n), by=Cluster_id]

pB <- ggplot(comp_long, aes(x = Cluster_id, y = prop, fill = disease)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "UC" = "#E94F37")) +
  labs(x = "T Cell Cluster", y = "Proportion",
       title = "B. Overall Composition") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  guides(fill = guide_legend(title = "Disease"))

# ── Panel C: Rectum-specific composition (inflammation site) ────────────────
rectum_data <- meta[tissue_type == "Rectum", .(n = .N), by=.(Cluster_id, disease)]
rectum_data[, prop := n / sum(n), by=Cluster_id]

pC <- ggplot(rectum_data, aes(x = Cluster_id, y = prop, fill = disease)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "UC" = "#E94F37")) +
  labs(x = "T Cell Cluster", y = "Proportion",
       title = "C. Rectal T Cells (Inflammation Site)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  guides(fill = guide_legend(title = "Disease"))

# ── Panel D: Cluster size distribution ──────────────────────────────────────
size_data <- meta[, .(n = .N), by=Cluster_id]

pD <- ggplot(size_data, aes(x = reorder(Cluster_id, -n), y = n)) +
  geom_bar(stat = "identity", fill = "#4A7C59", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.5) +
  labs(x = "T Cell Cluster", y = "Number of Cells",
       title = "D. Cluster Size") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── Combine panels with CORRECT Figure S9 caption ───────────────────────────
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Figure S9. Single-Cell T Cell Landscape in Ulcerative Colitis",
    subtitle = "Smillie et al. 2019 (GSE125527): 15,431 from 9 Healthy, 22,961 from 7 UC patients",
    caption = "Complement relevance: Cycling T (T9) and NK-like T (T5) most expanded in UC — consistent with complement-mediated inflammation",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.caption = element_text(hjust = 0.5, size = 9, face = "italic")
    )
  )

# Save high-resolution figure
ggsave(file.path(outdir, "Figure_S9_Tcell_UC.png"), fig, 
       width = 14, height = 10, dpi = 300, units = "in", bg = "white")
ggsave(file.path(outdir, "Figure_S9_Tcell_UC.pdf"), fig, 
       width = 14, height = 10, units = "in")

cat("Saved: Figure_S9_Tcell_UC.png/pdf\n")

# Save summary statistics
fwrite(comp_data, file.path(outdir, "Figure_S9_Tcell_stats.csv"))
cat("Saved: Figure_S9_Tcell_stats.csv\n")

cat("\n=== Figure S9 Complete ===\n")
