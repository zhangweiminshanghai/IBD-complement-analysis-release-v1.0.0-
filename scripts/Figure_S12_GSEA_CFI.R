# ============================================================================
# Figure S12: GSEA (GO Biological Process): CFI-high vs CFI-low CD samples
# Data Source: GSE16879 (bulk RNA-seq)
# Purpose: Re-generate compliant figure with correct "Figure S12" caption
# ============================================================================

library(data.table)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)

setwd("D:/IBD_Project")
outdir <- "D:/IBD_Project/final_suppl_figs"

# ── Load workspace with sample group information ────────────────────────────
if (file.exists("workspace_auto.RData")) {
  cat("Loading workspace with sample groups...\n")
  load("workspace_auto.RData")
} else {
  stop("workspace_auto.RData not found. Please run preprocessing first.")
}

# ── Get GSE16879 expression and sample groups ───────────────────────────────
# expr_16879_clean: 22168 genes x 133 samples
# group_train: named vector with CD/Control labels

cat("GSE16879 expression matrix:", nrow(expr_16879_clean), "genes x", 
    ncol(expr_16879_clean), "samples\n")

# Extract CD samples using pdata_16879_clean (has GSM IDs and disease status)
if (!exists("pdata_16879_clean")) {
  stop("pdata_16879_clean not found in workspace")
}

# Get GSM IDs from pdata rownames
pdata_gsm <- rownames(pdata_16879_clean)
cat("pdata GSM samples:", length(pdata_gsm), "\n")

# Find common samples between pdata and expression matrix
common_gsm <- intersect(pdata_gsm, colnames(expr_16879_clean))
cat("Common samples between pdata and expression:", length(common_gsm), "\n")

# Get disease status from pdata
disease_status <- pdata_16879_clean$disease[match(common_gsm, pdata_gsm)]
names(disease_status) <- common_gsm

cd_samples <- common_gsm[disease_status == "CD"]
control_samples <- common_gsm[disease_status == "Control"]

cat("CD samples:", length(cd_samples), "\n")
cat("Control samples:", length(control_samples), "\n")

# Subset expression matrix to CD samples
expr_cd <- expr_16879_clean[, cd_samples, drop=FALSE]

# ── Define CFI-high vs CFI-low groups in CD samples ─────────────────────────
cfi_expr <- as.numeric(expr_cd["CFI", ])
names(cfi_expr) <- colnames(expr_cd)

# Split by median CFI expression
cfi_median <- median(cfi_expr)
cfi_high <- names(cfi_expr)[cfi_expr > cfi_median]
cfi_low <- names(cfi_expr)[cfi_expr <= cfi_median]

cat("CFI-high CD samples:", length(cfi_high), "\n")
cat("CFI-low CD samples:", length(cfi_low), "\n")

# ── Differential expression analysis ────────────────────────────────────────
cat("Running differential expression...\n")
de_results <- data.frame(
  gene = rownames(expr_cd),
  logFC = NA,
  pvalue = NA
)

for (i in 1:nrow(expr_cd)) {
  high_expr <- as.numeric(expr_cd[i, cfi_high])
  low_expr <- as.numeric(expr_cd[i, cfi_low])
  
  # Log2 fold change
  de_results$logFC[i] <- mean(high_expr) - mean(low_expr)
  
  # T-test
  tt <- t.test(high_expr, low_expr)
  de_results$pvalue[i] <- tt$p.value
}

# Adjust p-values
de_results$padj <- p.adjust(de_results$pvalue, method = "fdr")

# Sort by logFC
de_results <- de_results[order(-de_results$logFC), ]

# Save DE results
write.csv(de_results, file.path(outdir, "Figure_S12_DE_CFI_high_vs_low.csv"), row.names = FALSE)
cat("Saved: Figure_S12_DE_CFI_high_vs_low.csv\n")

# ── Prepare gene list for GSEA ──────────────────────────────────────────────
geneList <- de_results$logFC
names(geneList) <- de_results$gene
geneList <- sort(geneList, decreasing = TRUE)

# Convert to ENTREZ IDs
gene_df <- bitr(de_results$gene, fromType = "SYMBOL", 
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Map ENTREZ IDs back to ranked list
entrez_list <- geneList[names(geneList) %in% gene_df$SYMBOL]
names(entrez_list) <- gene_df$ENTREZID[match(names(entrez_list), gene_df$SYMBOL)]
entrez_list <- entrez_list[!is.na(names(entrez_list))]
entrez_list <- sort(entrez_list, decreasing = TRUE)

# ── Run GSEA for GO Biological Process ──────────────────────────────────────
cat("Running GSEA for GO BP...\n")
gsea_go <- gseGO(
  geneList = entrez_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

cat("Significant GO BP terms:", nrow(as.data.frame(gsea_go)), "\n")

# Save GSEA results
if (nrow(as.data.frame(gsea_go)) > 0) {
  write.csv(as.data.frame(gsea_go), 
            file.path(outdir, "Figure_S12_GSEA_results.csv"), 
            row.names = FALSE)
  cat("Saved: Figure_S12_GSEA_results.csv\n")
}

# ── Create plots ────────────────────────────────────────────────────────────
# Panel A: GSEA dot plot
if (nrow(as.data.frame(gsea_go)) > 0) {
  pA <- dotplot(gsea_go, showCategory = 15) +
    labs(title = "A. GSEA Dot Plot (GO BP)") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
} else {
  pA <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "No significant GO terms") +
    theme_void() + labs(title = "A. GSEA Dot Plot (GO BP)")
}

# Panel B: GSEA ridge plot
if (nrow(as.data.frame(gsea_go)) > 0) {
  pB <- ridgeplot(gsea_go, showCategory = 15) +
    labs(title = "B. Enrichment Distribution") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))
} else {
  pB <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "No significant GO terms") +
    theme_void() + labs(title = "B. Enrichment Distribution")
}

# Panel C: Top pathways bar plot
if (nrow(as.data.frame(gsea_go)) > 0) {
  go_df <- as.data.frame(gsea_go)
  # Sort by NES (Normalized Enrichment Score)
  go_df <- go_df[order(-abs(go_df$NES)), ]
  top_go <- head(go_df, 15)
  top_go$Description <- factor(top_go$Description, levels = rev(top_go$Description))
  
  pC <- ggplot(top_go, aes(x = NES, y = Description)) +
    geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
    geom_text(aes(label = sprintf("p=%.2e", p.adjust)), hjust = -0.1, size = 2.5) +
    labs(x = "Normalized Enrichment Score (NES)", y = "GO Biological Process",
         title = "C. Top Enriched Pathways") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 9)) +
    xlim(min(top_go$NES) * 1.3, max(top_go$NES) * 1.3)
} else {
  pC <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "No significant GO terms") +
    theme_void() + labs(title = "C. Top Enriched Pathways")
}

# Panel D: GSEA running score (top pathway)
if (nrow(as.data.frame(gsea_go)) > 0) {
  top_pathway <- go_df$ID[1]
  pD <- gseaplot2(gsea_go, geneSetID = top_pathway, 
                  title = paste("D. Running Score:", go_df$Description[1])) +
    theme_bw(base_size = 11)
} else {
  pD <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "No significant GO terms") +
    theme_void() + labs(title = "D. Running Score")
}

# ── Combine panels with CORRECT Figure S12 caption ──────────────────────────
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Figure S12. GSEA (GO Biological Process): CFI-high vs CFI-low CD samples (GSE16879)",
    subtitle = sprintf("Complement-mediated inflammation pathways enriched in CFI-high Crohn's disease (n=%d CFI-high, n=%d CFI-low)",
                       length(cfi_high), length(cfi_low)),
    caption = "GSEA performed on CD samples split by median CFI expression; NES = Normalized Enrichment Score",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.caption = element_text(hjust = 0.5, size = 9, face = "italic")
    )
  )

# Save high-resolution figure
ggsave(file.path(outdir, "Figure_S12_GSEA_CFI.png"), fig, 
       width = 16, height = 12, dpi = 300, units = "in", bg = "white")
ggsave(file.path(outdir, "Figure_S12_GSEA_CFI.pdf"), fig, 
       width = 16, height = 12, units = "in")

cat("Saved: Figure_S12_GSEA_CFI.png/pdf\n")
cat("\n=== Figure S12 Complete ===\n")
