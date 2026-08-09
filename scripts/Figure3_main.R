###############################################################################
# Figure 3 (MAIN) - CFI associates with myeloid immune infiltration and a
# dual innate-adaptive co-infiltration architecture.
#
# Panels:
#   (A) CIBERSORT violin plot of 22 immune cell fractions, CD vs Control
#   (B) CFI expression vs immune fractions (boxplots by Low/High fraction,
#       Spearman rho annotated per fraction)
#   (C) Spearman correlation heatmap of the 22 immune cell fractions
#   (D) Innate co-infiltration module (Neutrophils, Mast cells activated, Monocytes)
#   (E) Adaptive/regulatory module (Tregs, CD4 memory resting, T follicular helper)
#
# Cohorts combined: GSE16879 (85: 73 CD/12 Control) + GSE75214 (97: 75 CD/22 Control)
###############################################################################

suppressMessages({
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(tidyr)
  library(dplyr)
  library(scales)
  library(grid)
})

setwd("D:/IBD_Project")
dir.create("final_main_figs", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Load data
## ---------------------------------------------------------------------------
load("workspace_auto.RData")   # provides pdata_16879_clean, expr_16879_clean, expr_75214_clean

frac_cols <- c("B cells naive","B cells memory","Plasma cells","T cells CD8",
               "T cells CD4 naive","T cells CD4 memory resting",
               "T cells CD4 memory activated","T cells follicular helper",
               "T cells regulatory (Tregs)","T cells gamma delta",
               "NK cells resting","NK cells activated","Monocytes",
               "Macrophages M0","Macrophages M1","Macrophages M2",
               "Dendritic cells resting","Dendritic cells activated",
               "Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")

cib1 <- read.csv("CIBERSORT_Results_GSE16879.csv", row.names = 1, check.names = FALSE)[, frac_cols]
cib2 <- read.csv("CIBERSORT_Results_GSE75214.csv", row.names = 1, check.names = FALSE)[, frac_cols]

# Disease labels
lab1 <- setNames(as.character(pdata_16879_clean$group), rownames(pdata_16879_clean))  # GSE16879
sg   <- read.csv("A5_A8_MR/GSE75214_sample_groups.csv")
lab2 <- setNames(as.character(sg$group), sg$GSM)                                       # GSE75214

s1 <- intersect(rownames(cib1), names(lab1))
s2 <- intersect(rownames(cib2), names(lab2))
cib1 <- cib1[s1, ]; cib2 <- cib2[s2, ]

# CFI expression (per-cohort z-score to allow pooling across batches)
cfi1 <- as.numeric(expr_16879_clean["CFI", s1]); names(cfi1) <- s1
cfi2 <- as.numeric(expr_75214_clean["CFI", s2]); names(cfi2) <- s2
cfi1z <- scale(cfi1)[,1]; cfi2z <- scale(cfi2)[,1]

# Combined data frame
comb <- rbind(
  data.frame(Sample = s1, Cohort = "GSE16879", Group = lab1[s1],
             CFIz = cfi1z[s1], cib1, check.names = FALSE),
  data.frame(Sample = s2, Cohort = "GSE75214", Group = lab2[s2],
             CFIz = cfi2z[s2], cib2, check.names = FALSE)
)
comb$Group <- factor(comb$Group, levels = c("Control","CD"))
frac_mat   <- as.matrix(comb[, frac_cols])   # 182 x 22
cat(sprintf("Combined samples: %d (CD=%d, Control=%d)\n",
            nrow(comb), sum(comb$Group=="CD"), sum(comb$Group=="Control")))

## ---------------------------------------------------------------------------
## Aesthetics
## ---------------------------------------------------------------------------
grp_cols <- c(Control = "#4C72B0", CD = "#C44E52")
theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 14),
        legend.key = element_blank())

## ---------------------------------------------------------------------------
## Panel A - Violin of 22 fractions, CD vs Control
## ---------------------------------------------------------------------------
dfA <- comb %>%
  pivot_longer(all_of(frac_cols), names_to = "CellType", values_to = "Fraction")
dfA$CellType <- factor(dfA$CellType, levels = frac_cols)

# Wilcoxon significance per cell type (CD vs Control)
sigA <- dfA %>% group_by(CellType) %>%
  summarise(p = tryCatch(wilcox.test(Fraction ~ Group)$p.value, error = function(e) NA_real_),
            ymax = max(Fraction, na.rm = TRUE), .groups = "drop") %>%
  mutate(star = cut(p, c(-Inf,0.001,0.01,0.05,Inf), labels = c("***","**","*","")))

pA <- ggplot(dfA, aes(CellType, Fraction, fill = Group)) +
  geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.55,
              position = position_dodge(0.8)) +
  geom_boxplot(width = 0.18, outlier.size = 0.25, alpha = 0.9,
               position = position_dodge(0.8), colour = "grey25", linewidth = 0.3) +
  geom_text(data = sigA, aes(x = CellType, y = ymax + 0.03, label = star),
            inherit.aes = FALSE, size = 4.2, colour = "grey20") +
  scale_fill_manual(values = grp_cols, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = NULL, y = "Estimated fraction",
       title = "(A) CIBERSORT immune cell fractions: CD vs Control (both cohorts)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "top")

## ---------------------------------------------------------------------------
## Panel B - CFI expression vs immune fractions (boxplots + Spearman rho)
## ---------------------------------------------------------------------------
# Spearman rho of CFI (z) vs each fraction across the combined cohort
rho_tab <- data.frame(CellType = frac_cols,
                      rho = sapply(frac_cols, function(cc)
                        suppressWarnings(cor(comb$CFIz, comb[[cc]], method = "spearman"))),
                      p   = sapply(frac_cols, function(cc)
                        suppressWarnings(cor.test(comb$CFIz, comb[[cc]], method = "spearman")$p.value)),
                      row.names = NULL)
rho_tab <- rho_tab[order(-abs(rho_tab$rho)), ]
write.csv(rho_tab, "final_main_figs/Figure3_CFI_fraction_Spearman.csv", row.names = FALSE)

# Low/High split per fraction (median), y = CFI z-score
dfB <- comb %>%
  pivot_longer(all_of(frac_cols), names_to = "CellType", values_to = "Fraction") %>%
  group_by(CellType) %>%
  mutate(Level = ifelse(Fraction > median(Fraction, na.rm = TRUE), "High", "Low")) %>%
  ungroup()
dfB$Level <- factor(dfB$Level, levels = c("Low","High"))

# Facet strip labels use short names + Spearman rho + stars, ordered by |rho|
rho_tab$star <- cut(rho_tab$p, c(-Inf,0.001,0.01,0.05,Inf), labels = c("***","**","*","ns"))
short_names <- c(
  `B cells naive`="B naive", `B cells memory`="B mem", `Plasma cells`="Plasma",
  `T cells CD8`="CD8 T", `T cells CD4 naive`="CD4 naive",
  `T cells CD4 memory resting`="CD4 mem rest", `T cells CD4 memory activated`="CD4 mem act",
  `T cells follicular helper`="Tfh", `T cells regulatory (Tregs)`="Tregs",
  `T cells gamma delta`="g/d T", `NK cells resting`="NK rest", `NK cells activated`="NK act",
  `Monocytes`="Mono", `Macrophages M0`="M0", `Macrophages M1`="M1", `Macrophages M2`="M2",
  `Dendritic cells resting`="DC rest", `Dendritic cells activated`="DC act",
  `Mast cells resting`="Mast rest", `Mast cells activated`="Mast act",
  `Eosinophils`="Eos", `Neutrophils`="Neutro")
lab_map <- setNames(sprintf("%s\n(\u03c1=%.2f%s)", short_names[rho_tab$CellType], rho_tab$rho,
                            ifelse(rho_tab$star=="ns","", as.character(rho_tab$star))),
                    rho_tab$CellType)
dfB$CellType <- factor(dfB$CellType, levels = rho_tab$CellType)

pB <- ggplot(dfB, aes(Level, CFIz, fill = Level)) +
  geom_boxplot(outlier.size = 0.2, width = 0.65, linewidth = 0.3, colour = "grey25") +
  facet_wrap(~ CellType, ncol = 11, labeller = labeller(CellType = lab_map)) +
  scale_fill_manual(values = c(Low = "#8DA0CB", High = "#E78AC3"), guide = "none") +
  labs(x = "Immune fraction level (median split)", y = "CFI expression (z-score)",
       title = "(B) CFI expression vs immune fractions (Spearman \u03c1 per fraction)") +
  theme_pub +
  theme(strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 9.2, lineheight = 0.9),
        axis.text = element_text(size = 8),
        panel.spacing = unit(0.35, "lines"))

## ---------------------------------------------------------------------------
## Panel C - Spearman correlation heatmap of 22 fractions
## ---------------------------------------------------------------------------
cmat <- cor(frac_mat, method = "spearman")
ord  <- hclust(as.dist(1 - cmat))$order
cmat <- cmat[ord, ord]
dfC  <- melt(cmat, varnames = c("V1","V2"), value.name = "rho")
dfC$V1 <- factor(dfC$V1, levels = rownames(cmat))
dfC$V2 <- factor(dfC$V2, levels = rownames(cmat))

pC <- ggplot(dfC, aes(V1, V2, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-1,1), name = "Spearman \u03c1") +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "(C) Immune cell fraction correlation (Spearman)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.3),
        axis.text.y = element_text(size = 8.3),
        panel.grid = element_blank(),
        legend.position = "right")

## ---------------------------------------------------------------------------
## Panels D & E - co-infiltration modules (sub-heatmaps with rho labels)
## ---------------------------------------------------------------------------
module_heat <- function(members, title, subtitle) {
  m  <- cor(frac_mat[, members], method = "spearman")
  d  <- melt(m, varnames = c("V1","V2"), value.name = "rho")
  d$V1 <- factor(d$V1, levels = members); d$V2 <- factor(d$V2, levels = rev(members))
  ggplot(d, aes(V1, V2, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3.6,
              colour = ifelse(abs(d$rho) > 0.6, "white", "grey15")) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                         midpoint = 0, limits = c(-1,1), name = "\u03c1") +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5),
          axis.text.y = element_text(size = 8.5),
          panel.grid = element_blank(),
          plot.subtitle = element_text(size = 9, colour = "grey35"),
          legend.position = "right")
}

innate_mods   <- c("Neutrophils","Mast cells activated","Monocytes")
adaptive_mods <- c("T cells regulatory (Tregs)","T cells CD4 memory resting",
                   "T cells follicular helper")

pD <- module_heat(innate_mods,  "(D) Innate co-infiltration module",
                   "Neutrophils / Activated mast cells / Monocytes")
pE <- module_heat(adaptive_mods, "(E) Adaptive-regulatory module",
                   "Tregs / CD4 memory resting / T follicular helper")

## ---------------------------------------------------------------------------
## Compose
## ---------------------------------------------------------------------------
design <- "
AAAAAA
BBBBBB
CCCCDD
CCCCEE
"
fig <- pA + pB + pC + pD + pE +
  plot_layout(design = design, heights = c(1.05, 0.95, 1.05, 1.05))

title_txt <- paste0("Figure 3. CFI associates with myeloid immune infiltration and a ",
                    "dual innate\u2013adaptive co-infiltration architecture.")

fig <- fig + plot_annotation(
  title = title_txt,
  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0,
                                          margin = margin(b = 6)),
                plot.background = element_rect(fill = "white", colour = NA)))

## ---------------------------------------------------------------------------
## Save (300 DPI)
## ---------------------------------------------------------------------------
png_path <- "final_main_figs/Figure3_main.png"
pdf_path <- "final_main_figs/Figure3_main.pdf"

ggsave(png_path, fig, width = 16, height = 12, units = "in", dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave(pdf_path, fig, width = 16, height = 12, units = "in",
       bg = "white", limitsize = FALSE, device = "pdf")

cat("\n=== Top CFI-correlated fractions (Spearman) ===\n")
print(head(rho_tab[, c("CellType","rho","p","star")], 8), row.names = FALSE)
cat("\nSaved:\n ", normalizePath(png_path), "\n ", normalizePath(pdf_path), "\n")
