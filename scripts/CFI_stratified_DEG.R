#!/usr/bin/env Rscript
# ============================================================================
# CFI_stratified_DEG.R
# ----------------------------------------------------------------------------
# Purpose : Expression-stratified (median-split) differential expression to
#           support the "CFI as a downstream compensatory node" hypothesis
#           (Figure 5 / CFI-pathway analysis in the manuscript).
# Inputs  : data/processed/GSE16879_expr_normalized.txt  (genes x samples, log2)
#           data/processed/GSE16879_phenotype.csv         (sample, group)
# Outputs : results/CFI_high_vs_low_DEG.csv   (limma top-table)
#           results/figures/CFI_pathway_volcano.png/.pdf
#           results/CFI_pathway_correlations.csv
# Key params : split = median(CFI expression); |logFC|>0.5 & adj.P.Val<0.05
# Approx runtime : < 1 min (discovery cohort only)
# Required : limma, statmod, ggplot2
# ============================================================================

set.seed(2024)
suppressMessages({
  library(limma)
  library(ggplot2)
})

root   <- ".."
expr_f <- file.path(root, "data", "processed", "GSE16879_expr_normalized.txt")
phen_f <- file.path(root, "data", "processed", "GSE16879_phenotype.csv")
out_dir <- file.path(root, "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(expr_f)) {
  message("ERROR: processed matrix not found at ", expr_f,
          ". Place GSE16879 normalized expression there (see data/README.md).")
  quit(status = 1)
}

expr  <- as.matrix(read.delim(expr_f, check.names = FALSE, stringsAsFactors = FALSE))
pheno <- read.csv(phen_f, stringsAsFactors = FALSE)
# keep CD samples only (CFI stratification is meaningful in disease context)
pheno <- pheno[pheno$group %in% c("CD", "Crohn"), , drop = FALSE]
common <- intersect(colnames(expr), pheno$sample)
expr <- expr[, common, drop = FALSE]
pheno <- pheno[match(common, pheno$sample), ]

cfi <- expr["CFI", , drop = FALSE]
if (is.null(cfi) || length(cfi) == 0) {
  message("ERROR: CFI not present in expression matrix.")
  quit(status = 1)
}
cfi_vec <- as.numeric(cfi)

# Median split into CFI-high / CFI-low
med <- median(cfi_vec)
grp <- ifelse(cfi_vec >= med, "CFI_high", "CFI_low")
message(sprintf("CFI median = %.3f; high n=%d, low n=%d",
                med, sum(grp == "CFI_high"), sum(grp == "CFI_low")))

design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(factor(grp))
contr <- makeContrasts(CFI_high - CFI_low, levels = design)

fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit, trend = TRUE)
tt  <- topTable(fit, n = Inf, sort.by = "P")

tt$Change <- ifelse(tt$adj.P.Val < 0.05 & tt$logFC > 0.5, "up",
              ifelse(tt$adj.P.Val < 0.05 & tt$logFC < -0.5, "down", "ns"))
write.csv(tt, file.path(out_dir, "CFI_high_vs_low_DEG.csv"))

# CFI-pathway correlation: correlate CFI expression with the top DEGs / CFI pathway genes
candi <- rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5]
candi <- setdiff(candi, "CFI")
cor_out <- data.frame(gene = character(), pearson_r = numeric(),
                      p_value = numeric(), stringsAsFactors = FALSE)
for (g in candi[1:min(50, length(candi))]) {
  sub <- expr[g, , drop = FALSE]
  r <- cor.test(as.numeric(sub), cfi_vec, method = "pearson")
  cor_out <- rbind(cor_out, data.frame(gene = g,
                                        pearson_r = round(r$estimate, 4),
                                        p_value = format.pval(r$p.value, digits = 3),
                                        stringsAsFactors = FALSE))
}
write.csv(cor_out, file.path(out_dir, "CFI_pathway_correlations.csv"), row.names = FALSE)

# Volcano
p <- ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val), colour = Change)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_colour_manual(values = c(up = "#d73027", down = "#1a9850", ns = "#999999")) +
  theme_bw(base_size = 9) +
  labs(title = "CFI-high vs CFI-low (median split, GSE16879 CD)",
       x = "log2 FC", y = "-log10 adj.P") +
  geom_hline(yintercept = -log10(0.05), linetype = 2)
ggsave(file.path(fig_dir, "CFI_pathway_volcano.png"), p, width = 5, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "CFI_pathway_volcano.pdf"), p, width = 5, height = 4)

message("CFI_stratified_DEG.R complete: ",
        sum(tt$Change != "ns"), " significant DEGs (|logFC|>0.5, adj.P<0.05).")
print(sessionInfo())
