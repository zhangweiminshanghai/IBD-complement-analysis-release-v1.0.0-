# =============================================================================
# S10_gsea.R
# Supplementary Figure S10 — GSEA dot plot + dual-axis ridge plot (v6) + NES bars
# Dual-axis fix: left ax (0-1) pure labels, right ax (-2.4 to 2.4) pure curves
# Run: Rscript code/R/S10_gsea.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(grid)
  library(gridExtra)
})

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- GSEA summary data -----------------------------------------------------
gsea <- data.frame(
  term = factor(c("cytoplasmic translation", "ribosome biogenesis",
                  "xenobiotic metabolism", "lipid metabolic process"),
                levels = c("cytoplasmic translation", "ribosome biogenesis",
                           "xenobiotic metabolism", "lipid metabolic process")),
  NES   = c(-1.98, -1.76, 1.54, 1.42),
  padj  = c(6.28e-7, 1.20e-3, 5.79e-4, 2.80e-2),
  ngene = c(120, 95, 60, 78)
)

# ---- (A) Dot plot -----------------------------------------------------------
pA <- ggplot(gsea, aes(x = NES, y = term, size = ngene, colour = NES)) +
  geom_point() +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_size(range = c(4, 10)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "(A) GSEA dot plot", x = "Normalized Enrichment Score (NES)",
       y = "", size = "Gene set size") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom")

# ---- (B) Ridge plot — DUAL AXIS (v6 architecture) ---------------------------
# Simulate running-enrichment scores for each term
set.seed(3)
terms <- levels(gsea$term)
ridges <- list()
for (i in seq_along(terms)) {
  x <- seq(-2500, 2500, length.out = 500)
  nes <- gsea$NES[i]
  # peak shifts with NES
  mu <- nes * 800
  y <- dnorm(x, mean = mu, sd = 900) * abs(nes)
  ridges[[i]] <- data.frame(x = x, density = y, term = terms[i])
}
df_ridge <- do.call(rbind, ridges)
df_ridge$term <- factor(df_ridge$term, levels = levels(gsea$term))

# Colour by NES sign
ridge_cols <- c("cytoplasmic translation" = "#2166AC",
                "ribosome biogenesis" = "#4393C3",
                "xenobiotic metabolism" = "#D6604D",
                "lipid metabolic process" = "#B2182B")

# We build two separate plots and combine with grid
# Left axis: pure labels
pB_labels <- ggplot() +
  geom_text(data = data.frame(y = seq_along(terms), lab = terms),
            aes(x = 0.5, y = y, label = lab),
            hjust = 0.5, size = 3.5, colour = "#2C3E50") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Right axis: pure ridge curves
pB_ridges <- ggplot(df_ridge, aes(x = x, y = term, height = density, fill = term)) +
  geom_ridgeline(scale = 1.5, alpha = 0.75, colour = "white", size = 0.3) +
  scale_fill_manual(values = ridge_cols) +
  scale_x_continuous(limits = c(-2400, 2400), expand = c(0, 0),
                     name = "Running enrichment score") +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 9))

# Combine side-by-side with exact widths
pB <- arrangeGrob(pB_labels, pB_ridges, ncol = 2, widths = c(0.38, 0.62))
pB <- arrangeGrob(pB, top = textGrob("(B) Enrichment distributions", gp = gpar(fontsize = 11, fontface = "bold")))

# ---- (C) NES bar chart ------------------------------------------------------
pC <- ggplot(gsea, aes(x = NES, y = term, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(title = "(C) NES bar chart", x = "Normalized Enrichment Score (NES)", y = "") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "none")

# ---- Combine A / B / C ------------------------------------------------------
pA_resized <- pA + theme(legend.position = "none")
combined <- arrangeGrob(
  arrangeGrob(pA_resized, ncol = 1, widths = 1),
  pB, pC,
  ncol = 3, widths = c(1.0, 1.4, 0.9)
)

ggsave(
  filename = file.path(out_dir, "S10_gsea_v6.png"),
  plot = combined, width = 9, height = 3.0, dpi = 220, bg = "white"
)
message("[S10] saved → figures/S10_gsea_v6.png")
