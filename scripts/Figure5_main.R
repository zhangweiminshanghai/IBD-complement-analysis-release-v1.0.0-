# Figure 5: Expression-stratified transcriptomic analysis of CFI
# Publication-quality multi-panel figure (A-D) for Crohn's disease manuscript
# Panels: A=Volcano; B=GO enrichment (downregulated); C=GO enrichment (upregulated); D=Top DEG heatmap
# NOTE: KEGG pathway analysis requires the KEGG.db package + live KEGG REST access,
#       which are unavailable in this offline environment; GO (Biological Process)
#       enrichment is used for both pathway panels so the figure matches its caption.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(gridExtra)
  library(png)
})

output_dir <- "D:/IBD_Project/final_main_figs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================================
# PANEL A: Volcano plot of DEGs (CFI-high vs CFI-low)
# ============================================================================
cat("Loading DEG data...\n")
deg_data <- read.csv("D:/IBD_Project/final_suppl_figs/Figure_S12_DE_CFI_high_vs_low.csv",
                     stringsAsFactors = FALSE)

deg_data <- deg_data %>%
  mutate(log2FC = log2(abs(logFC)) * sign(logFC)) %>%
  filter(!is.na(log2FC) & !is.na(padj) & !is.infinite(log2FC))

padj_threshold <- 0.05
log2fc_threshold <- 1

deg_data <- deg_data %>%
  mutate(significance = case_when(
    padj < padj_threshold & log2FC > log2fc_threshold ~ "Up",
    padj < padj_threshold & log2FC < -log2fc_threshold ~ "Down",
    TRUE ~ "NS"))

n_up <- sum(deg_data$significance == "Up", na.rm = TRUE)
n_down <- sum(deg_data$significance == "Down", na.rm = TRUE)
n_total <- n_up + n_down
cat("DEG counts:", n_up, "up,", n_down, "down,", n_total, "total\n")

label_genes <- c("CFB", "C2", "REG1A", "REG1B", "GBP2", "GBP4")

volcano_plot <- ggplot(deg_data, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.5, size = 1.2) +
  scale_color_manual(
    values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
    labels = c("Up" = paste0("Up (n=", n_up, ")"),
               "Down" = paste0("Down (n=", n_down, ")"),
               "NS" = "Not significant"),
    name = NULL) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_threshold),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(data = subset(deg_data, gene %in% label_genes),
            aes(label = gene), size = 3, fontface = "bold", vjust = -0.5,
            check_overlap = TRUE, color = "black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        axis.title = element_text(size = 10), axis.text = element_text(size = 9),
        legend.position = "right", panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  labs(title = "A. Volcano Plot",
       x = expression(log[2](fold~change)),
       y = expression(-log[10](adjusted~italic(p)))) +
  coord_cartesian(xlim = c(min(deg_data$log2FC) * 1.05, max(deg_data$log2FC) * 1.05))

# ============================================================================
# PANEL B: GO enrichment of DOWNREGULATED genes
# ============================================================================
cat("\nPreparing Panel B: GO enrichment of downregulated genes...\n")
down_genes <- deg_data %>% filter(padj < 0.05, log2FC < -1) %>% pull(gene)
cat("Downregulated genes:", length(down_genes), "\n")
go_down_plot <- NULL
if (length(down_genes) > 0) {
  conv_down <- tryCatch(bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db), error = function(e) { cat("bitr err:", conditionMessage(e), "\n"); NULL })
  if (!is.null(conv_down) && nrow(conv_down) > 0) {
    ego_down <- tryCatch(enrichGO(gene = unique(conv_down$ENTREZID), OrgDb = org.Hs.eg.db,
                                  ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
                         error = function(e) { cat("enrichGO down err:", conditionMessage(e), "\n"); NULL })
    if (!is.null(ego_down) && nrow(ego_down@result) > 0) {
      go_down_plot <- dotplot(ego_down, showCategory = 15) +
        ggtitle("B. GO Enrichment (Downregulated Genes)") +
        theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
              axis.title = element_text(size = 10), axis.text = element_text(size = 9),
              axis.text.y = element_text(size = 7), legend.position = "right",
              panel.grid.major = element_line(color = "grey90"), panel.grid.minor = element_blank())
      cat("GO down terms:", nrow(ego_down@result), "\n")
    }
  }
}
if (is.null(go_down_plot)) {
  go_down_plot <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "GO enrichment\n(unavailable)", size = 4, hjust = 0.5) +
    theme_void() + labs(title = "B. GO Enrichment (Downregulated)") +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
}

# ============================================================================
# PANEL C: GO enrichment of UPREGULATED genes (pathway analysis)
# ============================================================================
cat("\nPreparing Panel C: GO enrichment of upregulated genes...\n")
up_genes <- deg_data %>% filter(padj < 0.05, log2FC > 1) %>% pull(gene)
cat("Upregulated genes:", length(up_genes), "\n")
go_up_plot <- NULL
if (length(up_genes) > 0) {
  conv_up <- tryCatch(bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db), error = function(e) { cat("bitr err:", conditionMessage(e), "\n"); NULL })
  if (!is.null(conv_up) && nrow(conv_up) > 0) {
    ego_up <- tryCatch(enrichGO(gene = unique(conv_up$ENTREZID), OrgDb = org.Hs.eg.db,
                                ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
                       error = function(e) { cat("enrichGO up err:", conditionMessage(e), "\n"); NULL })
    if (!is.null(ego_up) && nrow(ego_up@result) > 0) {
      go_up_plot <- dotplot(ego_up, showCategory = 15) +
        ggtitle("C. GO Enrichment (Upregulated Genes)") +
        theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
              axis.title = element_text(size = 10), axis.text = element_text(size = 9),
              axis.text.y = element_text(size = 7), legend.position = "right",
              panel.grid.major = element_line(color = "grey90"), panel.grid.minor = element_blank())
      cat("GO up terms:", nrow(ego_up@result), "\n")
    }
  }
}
if (is.null(go_up_plot)) {
  go_up_plot <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "GO enrichment\n(unavailable)", size = 4, hjust = 0.5) +
    theme_void() + labs(title = "C. GO Enrichment (Upregulated)") +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
}

# ============================================================================
# PANEL D: Heatmap of top DEGs
# ============================================================================
cat("\nPreparing heatmap data...\n")
expr_data <- read.table("D:/IBD_Project/GSE16879_expr_normalized.txt", header = TRUE,
                        sep = "\t", row.names = 1, stringsAsFactors = FALSE,
                        check.names = FALSE)

top_genes <- deg_data %>%
  filter(gene != "CFI") %>%
  arrange(padj) %>%
  filter(gene %in% rownames(expr_data)) %>%
  head(25) %>%
  pull(gene)

heatmap_matrix <- expr_data[rownames(expr_data) %in% top_genes, ]
top_genes_found <- rownames(heatmap_matrix)
cat("Top genes for heatmap:", length(top_genes_found), "\n")

if ("CFI" %in% rownames(expr_data)) {
  cfi_expr <- as.numeric(expr_data["CFI", ])
  names(cfi_expr) <- colnames(expr_data)
  cfi_median <- median(cfi_expr, na.rm = TRUE)
  sample_groups <- ifelse(cfi_expr >= cfi_median, "CFI-high", "CFI-low")
  names(sample_groups) <- colnames(expr_data)
}

common_samples <- intersect(colnames(heatmap_matrix), names(sample_groups))
heatmap_matrix <- heatmap_matrix[, common_samples]
sample_groups <- sample_groups[common_samples]
sample_order <- order(sample_groups)
heatmap_matrix <- heatmap_matrix[, sample_order]
sample_groups <- sample_groups[sample_order]

heatmap_zscore <- t(scale(t(heatmap_matrix)))
heatmap_zscore[heatmap_zscore > 3] <- 3
heatmap_zscore[heatmap_zscore < -3] <- -3

annotation_col <- data.frame(
  Group = factor(sample_groups, levels = c("CFI-low", "CFI-high")),
  row.names = colnames(heatmap_matrix))
ann_colors <- list(Group = c("CFI-high" = "#E41A1C", "CFI-low" = "#377EB8"))

heatmap_png <- file.path(output_dir, "Figure5_panelD_heatmap.png")
png(heatmap_png, width = 1000, height = 900, res = 100)
pheatmap(heatmap_zscore, cluster_rows = TRUE, cluster_cols = FALSE,
         show_colnames = FALSE, show_rownames = TRUE,
         annotation_col = annotation_col, annotation_colors = ann_colors,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         fontsize_row = 9, fontsize = 10, main = "D. Top DEGs Heatmap",
         border_color = NA, scale = "none")
dev.off()
cat("Heatmap created:", heatmap_png, "\n")

# ============================================================================
# COMPOSE FINAL FIGURE
# ============================================================================
cat("\nComposing final figure...\n")
volcano_png <- file.path(output_dir, "Figure5_panelA_volcano.png")
ggsave(volcano_png, volcano_plot, width = 8, height = 6, units = "in", dpi = 300)
goB_png <- file.path(output_dir, "Figure5_panelB_GO.png")
ggsave(goB_png, go_down_plot, width = 8, height = 6, units = "in", dpi = 300)
goC_png <- file.path(output_dir, "Figure5_panelC_GO.png")
ggsave(goC_png, go_up_plot, width = 8, height = 6, units = "in", dpi = 300)

img_a <- readPNG(volcano_png)
img_b <- readPNG(goB_png)
img_c <- readPNG(goC_png)
img_d <- readPNG(heatmap_png)

final_png <- file.path(output_dir, "Figure5_main.png")
png(final_png, width = 4800, height = 3600, res = 300)
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1, 0.5), oma = c(0, 0, 2, 0))
for (im in list(img_a, img_b, img_c, img_d)) {
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  rasterImage(im, 0, 0, 1, 1)
}
mtext("Figure 5. Expression-stratified transcriptomic analysis of CFI",
      outer = TRUE, cex = 1.2, font = 2, line = 0.5)
dev.off()
cat("Final PNG created:", final_png, "\n")

final_pdf <- file.path(output_dir, "Figure5_main.pdf")
pdf(final_pdf, width = 16, height = 12)
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1, 0.5), oma = c(0, 0, 2, 0))
for (im in list(img_a, img_b, img_c, img_d)) {
  plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
  rasterImage(im, 0, 0, 1, 1)
}
mtext("Figure 5. Expression-stratified transcriptomic analysis of CFI",
      outer = TRUE, cex = 1.2, font = 2, line = 0.5)
dev.off()
cat("Final PDF created:", final_pdf, "\n")

temp_files <- list.files(output_dir, pattern = "^Figure5_panel", full.names = TRUE)
file.remove(temp_files)
cat("Temporary files cleaned up.\n")

# ============================================================================
# VERIFICATION
# ============================================================================
cat("\n========================================\n")
cat("FIGURE 5 GENERATION COMPLETE\n")
cat("========================================\n\n")
if (file.exists(final_png)) {
  img_dim <- dim(readPNG(final_png))
  cat("PNG file: ", final_png, "\n")
  cat("  Dimensions: ", img_dim[2], "x", img_dim[1], "pixels\n")
  cat("  DPI: 300\n  Size: ", img_dim[2]/300, "x", img_dim[1]/300, "inches\n")
}
cat("\nSummary:\n")
cat("  Total DEGs (padj<0.05, |log2FC|>1): ", n_total, "\n")
cat("  Upregulated: ", n_up, "\n")
cat("  Downregulated: ", n_down, "\n")
cat("  Genes for heatmap: ", length(top_genes_found), "\n")
cat("========================================\n")
