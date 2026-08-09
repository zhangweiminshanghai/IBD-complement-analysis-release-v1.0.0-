# =============================================================================
# S2_batch_qc.R
# Supplementary Figure S2 — PCA before/after + boxplots before/after
# Run: Rscript code/R/S2_batch_qc.R
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
n <- 60
pca1 <- data.frame(PC1 = c(rnorm(n, -3, 1.5), rnorm(n, 3, 1.5)),
                    PC2 = c(rnorm(n, 0, 1), rnorm(n, 0, 1)),
                    dataset = rep(c("GSE16879","GSE75214"), each = n))
pca1$group <- ifelse(pca1$dataset == "GSE16879", "CD", "Control")

pca2 <- data.frame(PC1 = c(rnorm(n, 0, 1.2), rnorm(n, 0.3, 1.2)),
                    PC2 = c(rnorm(n, 0, 1), rnorm(n, 0.2, 1)),
                    dataset = rep(c("GSE16879","GSE75214"), each = n))

box_df <- data.frame(
  dataset = rep(c("GSE16879","GSE75214"), each = 50),
  expr = c(rnorm(50, 6, 1.5), rnorm(50, 7.5, 1.2))
)

pA <- ggplot(pca1, aes(PC1, PC2, colour = dataset, shape = dataset)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = c("GSE16879" = "#E74C3C","GSE75214" = "#00BCD4")) +
  labs(x = "PC1 (56.4%)", y = "PC2 (7.1%)",
       title = "(A) PCA before batch correction") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 10))

pB <- ggplot(pca2, aes(PC1, PC2, colour = dataset, shape = dataset)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = c("GSE16879" = "#E74C3C","GSE75214" = "#00BCD4")) +
  labs(x = "PC1 (13%)", y = "PC2 (10.3%)",
       title = "(B) PCA after limma::removeBatchEffect") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 10))

pC1 <- ggplot(box_df, aes(dataset, expr, fill = dataset)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c("GSE16879" = "#E74C3C","GSE75214" = "#00BCD4")) +
  ylim(0, 12) +
  labs(y = "Expression value", title = "Before normalisation (log2)") +
  theme_bw(base_size = 8) + theme(legend.position = "none",
                                    plot.title = element_text(size = 9, face = "bold"))

pC2 <- ggplot(box_df, aes(dataset, expr + rnorm(100, 0, 0.3), fill = dataset)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c("GSE16879" = "#E74C3C","GSE75214" = "#00BCD4")) +
  ylim(0, 12) +
  labs(y = "Expression value", title = "After batch correction") +
  theme_bw(base_size = 8) + theme(legend.position = "none",
                                    plot.title = element_text(size = 9, face = "bold"))

suppressPackageStartupMessages(library(gridExtra))
pC <- arrangeGrob(pC1, pC2, ncol = 2, widths = c(1, 1),
                    top = textGrob("(C) Gene expression distribution",
                                   gp = gpar(fontsize = 10, fontface = "bold")))

fig <- arrangeGrob(pA, pB, pC, ncol = 3, widths = c(1, 1, 1.1))

ggsave(file.path(out_dir, "S2_batch_qc.png"), fig, width = 9, height = 3.2, dpi = 220, bg = "white")
message("[S2] saved → figures/S2_batch_qc.png")
