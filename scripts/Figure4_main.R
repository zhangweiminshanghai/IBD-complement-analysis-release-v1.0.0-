#!/usr/bin/env Rscript
# ==============================================================================
# Figure 4 — Crohn's disease single-cell manuscript
# Single-cell resolution reveals CFI-enriched inflammatory stromal cells in
# ileal CD and a KCNE3-mediated communication hub.
#
# Panels:
#   (A)  Cell-type UMAP
#   (B)  CFI  feature plot
#   (C)  C2   feature plot
#   (D)  Cell-type composition barplot
#   (E)  CFI+ fraction by disease location  (canonical + high-confidence)
#   (F)  UMAP highlighting ileal-CD high-confidence CFI+ population
#
# Data : Martin2019 / GSE134809  ~153 k cells
# Save : D:\IBD_Project\final_main_figs\Figure4_main.{png,pdf}
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(dplyr)
  library(tidyr)
})

message("[" , Sys.time() , "] Starting Figure 4 generation ...")

# --------------------------------------------------------------------------
# 0.  Paths
# --------------------------------------------------------------------------
RDS      <- "D:/IBD_Project/scRNA_data/Martin2019_results/Martin2019_processed_fixed.rds"
OUT_DIR  <- "D:/IBD_Project/final_main_figs"
OUT_PNG  <- file.path(OUT_DIR, "Figure4_main.png")
OUT_PDF  <- file.path(OUT_DIR, "Figure4_main.pdf")
FIG_W    <- 4800    # px at 300 dpi -> 16 inches
FIG_H    <- 3600    # px at 300 dpi -> 12 inches

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1.  Load Seurat object
# --------------------------------------------------------------------------
message("[" , Sys.time() , "] Loading RDS (~973 MB) ...")
t0 <- Sys.time()
sce <- readRDS(RDS)
message("[" , Sys.time() , "] Loaded in ", round(difftime(Sys.time(), t0, units = "secs"), 1), " s")
message("Dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells")

if (!inherits(sce, "Seurat")) {
  message("[INFO] Converting to Seurat ...")
  sce <- as.Seurat(sce, counts = "counts", data = "data")
}

# --------------------------------------------------------------------------
# 2.  Build disease-location annotation from orig.ident
# --------------------------------------------------------------------------
# From GEO GSE134809 sample records (full format per sample):
#
#  GSM3972009-2030  (22 samples, 2019 primary cohort):
#    Odd GSM IDs:  "Ileal Involved"    -> ileal CD
#    Even GSM IDs: "Ileal Uninvolved"  -> Control
#    GSM3972030_209 = patient 209 (odd) = Ileal Involved -> ileal CD
#
#  GSM4761136-1144  (9 samples, 2020 supplemental):
#    All "XX; PBMC" -> PBMC
#
# In the Seurat object: orig.ident = e.g. "GSM3972009_69"

md <- sce@meta.data
# Extract GSM number (integer) from e.g. "GSM3972009_69" -> 3972009
md$gsm_num <- as.integer(gsub("^GSM([0-9]+)_[0-9]+$", "\\1", md$orig.ident, perl = TRUE))

na_gsm <- sum(is.na(md$gsm_num))
if (na_gsm > 0) {
  message("[WARN] ", na_gsm, " cells with unrecognised GSM pattern")
}

md$disease_location <- "Unknown"
# Primary ileal cohort: GSM3972009-3972030 (odd = ileal CD, even = Control)
is_ileal <- !is.na(md$gsm_num) & md$gsm_num >= 3972009 & md$gsm_num <= 3972030
md$disease_location[is_ileal & (md$gsm_num %% 2) == 1]  <- "ileal CD"
md$disease_location[is_ileal & (md$gsm_num %% 2) == 0]  <- "Control"
# 2020 PBMC supplemental
md$disease_location[!is.na(md$gsm_num) & md$gsm_num >= 4761136 & md$gsm_num <= 4761144] <- "PBMC"

sce@meta.data <- md

cat("\n=== Disease-location annotation ===\n")
print(table(md$disease_location, useNA = "ifany"))
cat("\nPer group:\n")
print(md %>% group_by(disease_location) %>%
        summarise(n_cells = n(), n_samples = n_distinct(orig.ident)) %>% as.data.frame())

# --------------------------------------------------------------------------
# 3.  Cell type identities
# --------------------------------------------------------------------------
cat("\n=== Cell-type identities ===\n")
idents <- Seurat::Idents(sce)
cat("N levels: ", length(levels(idents)), "\n")

# --------------------------------------------------------------------------
# 4.  Gene presence check
# --------------------------------------------------------------------------
cat("\n=== Gene presence ===\n")
for (g in c("CFI","C2")) {
  cat("  ", g, if(g %in% rownames(sce)) "FOUND" else "NOT FOUND", "\n")
}

# --------------------------------------------------------------------------
# 5.  Compute CFI+ fractions
# --------------------------------------------------------------------------
cat("\n=== Computing CFI+ fractions ===\n")

tissue_groups <- c("ileal CD", "Control")

# CRITICAL: the @counts slot matrix is unnamed, so use positional indexing.
# cfi_expr is a plain numeric vector aligned with colnames(sce).
# md index cells by rownames (identical to colnames(sce)).
DefaultAssay(sce) <- "RNA"
rna_assay <- sce@assays$RNA
cfi_expr <- as.numeric(rna_assay@counts["CFI", ])
# Verify alignment with colnames
stopifnot(length(cfi_expr) == ncol(sce))
names(cfi_expr) <- colnames(sce)   # attach names so [named] indexing works

cat("CFI counts stats: min=", min(cfi_expr), " max=", max(cfi_expr),
    " median=", median(cfi_expr), " CFI+ (counts>0)=", sum(cfi_expr > 0), "\n")

# Top-decile threshold: among CFI+ cells (UMI > 0)
pos_cfi           <- cfi_expr[cfi_expr > 0]
top_decile_thresh <- quantile(pos_cfi, probs = 0.90, na.rm = TRUE)
cat("Top-decile threshold (high-confidence):", round(top_decile_thresh, 4), "\n")

results_df <- data.frame(
  Group         = character(),
  N_cells       = integer(),
  CFI_canonical = numeric(),
  CFI_highconf  = numeric(),
  stringsAsFactors = FALSE
)

for (grp in tissue_groups) {
  cells   <- rownames(md)[md$disease_location == grp]
  n_tot   <- length(cells)
  if (n_tot == 0) { cat("  [WARN] No cells for:", grp, "\n"); next }

  cfi_sub    <- cfi_expr[cells]
  n_canon    <- sum(cfi_sub > 0, na.rm = TRUE)
  n_hiconf   <- sum(cfi_sub >= top_decile_thresh, na.rm = TRUE)
  pct_canon  <- 100 * n_canon  / n_tot
  pct_hiconf <- 100 * n_hiconf / n_tot

  cat(sprintf("  %-12s  n=%7d  canonical=%6.2f%%  highconf=%6.2f%%\n",
              grp, n_tot, pct_canon, pct_hiconf))

  results_df <- rbind(results_df, data.frame(
    Group         = grp,
    N_cells       = n_tot,
    CFI_canonical = pct_canon,
    CFI_highconf  = pct_hiconf,
    stringsAsFactors = FALSE
  ))
}

cat("\n=== Manuscript reference ===\n")
cat("  ileal CD  canonical=6.1%  high-conf=0.89%\n")
cat("  Control   canonical=8.3%   high-conf=0.31%\n")
cat("\n=== Computed fractions ===\n")
print(results_df)

# --------------------------------------------------------------------------
# 6.  Colour palette for cell types
# --------------------------------------------------------------------------
idents <- Seurat::Idents(sce)
n_ct   <- length(levels(idents))
ct_pal <- setNames(
  hcl(seq(0, 360, length.out = n_ct + 1)[seq_len(n_ct)], l = 65, c = 100),
  levels(idents)
)

# --------------------------------------------------------------------------
# 7.  Build individual panels
# --------------------------------------------------------------------------

# ---- Panel A: Cell-type UMAP ----
message("\n[" , Sys.time() , "] Rendering Panel A ...")
pA <- DimPlot(sce, reduction = "umap", raster = TRUE, pt.size = 0.15,
              cols = ct_pal, label = TRUE, label.size = 3.5, repel = TRUE) +
  labs(title = "(A) Cell Type Identity") +
  theme_void(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", hjust = 0, margin = margin(b = 5)),
        legend.position  = "right", legend.justification = "top",
        legend.key.size   = unit(0.35, "lines"),
        legend.text       = element_text(size = 7))

# ---- Panel B: CFI feature plot ----
message("[" , Sys.time() , "] Rendering Panel B ...")
DefaultAssay(sce) <- "RNA"
pB <- FeaturePlot(sce, features = "CFI", reduction = "umap", raster = TRUE,
                  pt.size = 0.15, max.cutoff = "q95") +
  labs(title = "(B) CFI Expression") +
  theme_void(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", hjust = 0, margin = margin(b = 5)),
        legend.position  = "right", legend.justification = "top",
        legend.key.size   = unit(0.35, "lines"),
        legend.text       = element_text(size = 7))

# ---- Panel C: C2 feature plot ----
message("[" , Sys.time() , "] Rendering Panel C ...")
pC <- FeaturePlot(sce, features = "C2", reduction = "umap", raster = TRUE,
                  pt.size = 0.15, max.cutoff = "q95") +
  labs(title = "(C) C2 Expression") +
  theme_void(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", hjust = 0, margin = margin(b = 5)),
        legend.position  = "right", legend.justification = "top",
        legend.key.size   = unit(0.35, "lines"),
        legend.text       = element_text(size = 7))

# ---- Panel D: Cell-type composition (stacked bar by disease location) ----
message("[" , Sys.time() , "] Rendering Panel D ...")
tissue_md <- md[md$disease_location %in% tissue_groups, , drop = FALSE]
tissue_md$cell_type <- as.character(Seurat::Idents(sce)[rownames(tissue_md)])

comp_df <- tissue_md %>%
  group_by(disease_location, cell_type) %>%
  tally(name = "n") %>%
  group_by(disease_location) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()
colnames(comp_df)[1] <- "Location"

pD <- ggplot(comp_df, aes(x = Location, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65) +
  scale_fill_manual(values = ct_pal, name = "Cell Type") +
  labs(x = NULL, y = "Proportion (%)", title = "(D) Cell Composition") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x        = element_text(angle = 20, hjust = 1, size = 9),
        axis.text.y         = element_text(size = 9),
        plot.title          = element_text(face = "bold", hjust = 0),
        legend.key.size      = unit(0.3, "lines"),
        legend.text          = element_text(size = 6.5),
        panel.grid.minor     = element_blank(),
        panel.grid.major.x   = element_blank())

# ---- Panel E: CFI+ fraction barplot ----
message("[" , Sys.time() , "] Rendering Panel E ...")
results_df$Group <- factor(results_df$Group, levels = tissue_groups)
results_long <- pivot_longer(results_df,
                             cols = c(CFI_canonical, CFI_highconf),
                             names_to = "Threshold",
                             values_to = "Fraction")
results_long$Threshold <- ifelse(
  results_long$Threshold == "CFI_canonical",
  "Canonical (UMI > 0)", "High-conf. (top decile)"
)
results_long$Threshold <- factor(results_long$Threshold,
                                 levels = c("Canonical (UMI > 0)",
                                            "High-conf. (top decile)"))

pE <- ggplot(results_long, aes(x = Group, y = Fraction, fill = Threshold)) +
  geom_bar(position = position_dodge(preserve = "single"),
           stat = "identity", width = 0.60) +
  geom_text(aes(label = sprintf("%.2f%%", Fraction)),
            position = position_dodge(width = 0.60),
            vjust = -0.4, size = 2.8, color = "black") +
  scale_fill_manual(name = "Threshold",
    values = c("Canonical (UMI > 0)"    = "#2196F3",
               "High-conf. (top decile)" = "#E53935")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Disease Location", y = "CFI+ Cells (%)",
       title = "(E) CFI+ Fraction by Location") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x         = element_text(size = 9),
        axis.text.y          = element_text(size = 9),
        plot.title           = element_text(face = "bold", hjust = 0),
        legend.position       = "right", legend.justification = "top",
        legend.key.size       = unit(0.30, "lines"),
        legend.text           = element_text(size = 7),
        panel.grid.minor      = element_blank())

# ---- Panel F: High-confidence CFI+ cells highlighted on UMAP ----
message("[" , Sys.time() , "] Rendering Panel F ...")
hc_thresh  <- top_decile_thresh
ileal_mask <- md$disease_location == "ileal CD"
hiconf_cfi <- cfi_expr >= hc_thresh & ileal_mask

cat("High-conf CFI+ ileal cells:", sum(hiconf_cfi), "\n")

other_cells <- colnames(sce)[!hiconf_cfi]
hc_cells    <- colnames(sce)[hiconf_cfi]
set.seed(42)
keep_other  <- sample(other_cells, size = min(8000, length(other_cells)), replace = FALSE)
keep_cells  <- c(hc_cells, keep_other)
plot_sce    <- subset(sce, cells = keep_cells)

hi_flag <- ifelse(colnames(plot_sce) %in% hc_cells,
                  "High-conf.\nCFI+ (ileal CD)", "Other")
hi_flag <- factor(hi_flag, levels = c("Other", "High-conf.\nCFI+ (ileal CD)"))
plot_sce@meta.data$hi_flag <- hi_flag

pF <- DimPlot(plot_sce, reduction = "umap", raster = TRUE, pt.size = 0.22,
              group.by = "hi_flag",
              order = c("High-conf.\nCFI+ (ileal CD)"),
              cols = c("Other" = "#CCCCCC",
                       "High-conf.\nCFI+ (ileal CD)" = "#D32F2F")) +
  guides(fill = guide_legend(override.aes = list(size = 2))) +
  labs(title = "(F) High-conf. CFI+ (ileal CD)") +
  theme_void(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", hjust = 0, margin = margin(b = 5)),
        legend.position  = "right", legend.justification = "top",
        legend.key.size   = unit(0.35, "lines"),
        legend.text       = element_text(size = 7))

# --------------------------------------------------------------------------
# 8.  Figure title
# --------------------------------------------------------------------------
fig_title <- paste0(
  "Figure 4. Single-cell resolution reveals CFI-enriched inflammatory\n",
  "stromal cells in ileal CD and a KCNE3-mediated communication hub."
)

# --------------------------------------------------------------------------
# 9.  Patchwork layout — 3 columns × 2 rows grid
# --------------------------------------------------------------------------
message("\n[" , Sys.time() , "] Assembling figure ...")
fig <- wrap_plots(A = pA, B = pB, C = pC,
                  D = pD, E = pE, F = pF,
                  ncol = 3, nrow = 2) +
  plot_annotation(title = fig_title,
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12.5, hjust = 0,
                                    margin = margin(b = 8, l = 0, t = 4)),
      plot.margin   = margin(10, 10, 10, 10, "pt")
    ))

# --------------------------------------------------------------------------
# 10. Save PNG  (300 dpi, 4800 x 3600 px)
# --------------------------------------------------------------------------
message("[" , Sys.time() , "] Saving PNG ...")
t1 <- Sys.time()
ragg::agg_png(OUT_PNG, width = FIG_W, height = FIG_H, res = 300)
print(fig)
dev.off()
message("[" , Sys.time() , "] PNG saved in ",
        round(difftime(Sys.time(), t1, units = "secs"), 1), " s")

# --------------------------------------------------------------------------
# 11. Save PDF (vector)
# --------------------------------------------------------------------------
message("[" , Sys.time() , "] Saving PDF ...")
t2 <- Sys.time()
pdf(OUT_PDF, width = 16, height = 12)
print(fig)
dev.off()
message("[" , Sys.time() , "] PDF saved in ",
        round(difftime(Sys.time(), t2, units = "secs"), 1), " s")

# --------------------------------------------------------------------------
# 12. Verification
# --------------------------------------------------------------------------
message("\n=== Verification ===")
if (file.exists(OUT_PNG)) {
  fi <- file.info(OUT_PNG)
  message("PNG  : ", OUT_PNG)
  message("      Size  : ", round(fi$size / 1024^2, 2), " MB")
  message("      Date  : ", fi$mtime)
  img <- tryCatch(png::readPNG(OUT_PNG), error = function(e) NULL)
  if (!is.null(img)) {
    message("      Dims  : ", ncol(img), " x ", nrow(img), " px")
    message("      DPI   : 300 (requested)")
  }
} else {
  message("ERROR: PNG not found!")
}

if (file.exists(OUT_PDF)) {
  fi2 <- file.info(OUT_PDF)
  message("PDF  : ", OUT_PDF)
  message("      Size  : ", round(fi2$size / 1024^2, 2), " MB")
}

message("\n=== Final CFI+ fractions ===")
print(results_df)

message("\n[", Sys.time(), "] DONE.")
