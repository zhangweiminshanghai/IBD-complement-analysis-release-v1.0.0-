# =============================================================================
# S1_flowchart.R
# Supplementary Figure S1 — Integrative multi-omics & MR workflow
# *Gut* style: vertical 6-tier flowchart, dark-red core-narrative box
# Dependencies: ggplot2, grid, gridExtra
# Run: Rscript code/R/S1_flowchart.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# ---- colour palette ----------------------------------------------------------
tier_colours <- c(
  L1 = "#D6EAF8",  # light blue   – data acquisition
  L2 = "#FADBD8",  # light red    – DE / WGCNA
  L3 = "#E8DAEF",  # light purple – MR
  L4 = "#D5F5E3",  # light green  – scRNA-seq
  L5 = "#FCF3CF",  # light yellow – signature / docking
  L6 = "#D1F2EB"   # light cyan   – integration
)
border_colours <- c(
  L1 = "#2980B9", L2 = "#C0392B", L3 = "#8E44AD",
  L4 = "#27AE60", L5 = "#D4AC0D", L6 = "#16A085"
)

# ---- helper: one tier = one grob ---------------------------------------------
make_tier <- function(label, items, fill, border, width = 10, height = 1.6) {
  txt <- paste0("• ", items, collapse = "\n")
  grid.roundrect(
    width = unit(width, "cm"), height = unit(height, "cm"),
    r = unit(0.3, "cm"), gp = gpar(fill = fill, col = border, lwd = 1.8)
  )
  # text on top of rect
  grid.text(
    label = txt, x = 0.5, y = 0.5,
    gp = gpar(fontsize = 10, col = "#2C3E50", fontfamily = "sans")
  )
}

# We'll compose with ggplot grobs for cleaner control
tier_df <- data.frame(
  tier = c("L1", "L2a", "L2b", "L3", "L4", "L5a", "L5b", "L6"),
  y    = c(8, 6.8, 5.6, 4.4, 3.2, 2.0, 0.8, -0.4)
)
tier_df$fill   <- tier_colours[c("L1","L2","L2","L3","L4","L5","L5","L6")]
tier_df$border <- border_colours[c("L1","L2","L2","L3","L4","L5","L5","L6")]
tier_df$label  <- list(
  "Data acquisition & preprocessing\n• GSE16879 (discovery, n=85)\n• GSE75214 (validation, n=97)\n• Normalisation & batch correction (limma)",
  "Differential expression\n• limma-voom\n• |log2FC|>1, FDR<0.05\n• Up/down genes → pathway analysis",
  "WGCNA co-expression\n• Soft-threshold β=14\n• Yellow module\n• CFI, S100A8, S100A9, PAQR5, KCNE3",
  "Mendelian randomisation\n• C2 (protective) → CD\n• IVW OR=0.43 (95% CI 0.30–0.61)\n• p=2×10⁻⁶",
  "Single-cell RNA-seq\n• GSE134809 (Seurat)\n• 160,981 cells, 12 cell types\n• Cell-type-specific complement",
  "3-gene signature\n• LASSO → 10-fold CV\n• CFI · PAQR5 · KCNE3\n• AUC=0.90 (disc.); AUC=0.98 (val.)",
  "Structural pharmacology\n• Nafamostat mesylate docking\n• C1r/C1s vs CFI selectivity\n• Therapeutic implication",
  "Integration & validation\n• Bulk + scRNA + WGCNA + MR + ML\n• Cross-cohort bootstrap\n• Unified model"
)

# ---- base plot (invisible) ---------------------------------------------------
p_base <- ggplot(tier_df, aes(x = 0, y = y)) +
  geom_blank() +
  theme_void() +
  xlim(-6, 6) +
  coord_fixed(ratio = 0.5)

# ---- draw tiers as annotation_custom grobs -----------------------------------
grobs <- list()
for (i in seq_len(nrow(tier_df))) {
  g <- grobTree(
    rectGrob(
      x = 0, y = 0, width = 0.95, height = 0.85,
      gp = gpar(fill = tier_df$fill[i], col = tier_df$border[i], lwd = 1.8),
      name = paste0("rect_", i)
    ),
    textGrob(
      label = tier_df$label[[i]], x = 0, y = 0,
      gp = gpar(fontsize = 9, col = "#2C3E50"),
      name = paste0("txt_", i)
    ),
    vp = NULL
  )
  grobs[[i]] <- annotation_custom(g, xmin = -5, xmax = 5, ymin = tier_df$y[i] - 0.6, ymax = tier_df$y[i] + 0.6)
}

p <- p_base + grobs +
  geom_segment(aes(x = 0, xend = 0, y = 8.6, yend = -1.1),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               colour = "#7F8C8D", lwd = 1.2) +
  annotate("text", x = -4.5, y = tier_df$y, label = paste0("L", seq_len(nrow(tier_df))),
           size = 4, fontface = "bold", colour = "#2C3E50") +
  ggtitle("(A) Multi-omics analytical pipeline") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

# ---- core narrative box (separate grob) --------------------------------------
core_box <- grobTree(
  rectGrob(
    x = 0, y = 0, width = 0.7, height = 0.4,
    gp = gpar(fill = "#C0392B", col = "#922B21", lwd = 2.5)
  ),
  textGrob(
    label = "Core narrative:\nC2 initiated, CFI compensatory",
    x = 0, y = 0,
    gp = gpar(fontsize = 10, col = "white", fontface = "bold")
  )
)

p_final <- p + annotation_custom(core_box, xmin = 2.5, xmax = 5.5, ymin = 4.0, ymax = 4.8)

ggsave(
  filename = file.path(out_dir, "S1_flowchart.png"),
  plot = p_final, width = 9, height = 7, dpi = 220, bg = "white"
)
message("[S1] saved → figures/S1_flowchart.png")
