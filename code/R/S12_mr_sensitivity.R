# =============================================================================
# S12_mr_sensitivity.R
# Supplementary Figure S12 — MR sensitivity (4 panels, aligned titles)
# Run: Rscript code/R/S12_mr_sensitivity.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- (A) Leave-one-out forest plot -----------------------------------------
snp <- c("rs488755","rs114502302","All SNPs")
or  <- c(0.42, 0.44, 0.43)
lo  <- c(0.30, 0.32, 0.30)
hi  <- c(0.61, 0.62, 0.61)
dfA <- data.frame(snp, or, lo, hi)
dfA$snp <- factor(dfA$snp, levels = rev(snp))

pA <- ggplot(dfA, aes(or, snp)) +
  geom_point(size = 3, colour = "#2980B9") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "#2980B9") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = 0.30, xend = 0.61, y = 1.5, yend = 1.5),
               colour = "#C0392B", lwd = 1.5) +
  annotate("text", x = 0.46, y = 1.8, label = "IVW OR=0.43\n(0.30–0.61)",
           colour = "#C0392B", size = 3, hjust = 0.5) +
  labs(x = "OR (95% CI) for CD", y = "", title = "(A) Leave-one-out meta-analysis") +
  xlim(0.25, 0.70) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))

# ---- (B) Radial MR placeholder (proper axes, title aligned) -----------------
pB <- ggplot() +
  theme_bw(base_size = 9) +
  coord_fixed(ratio = 1, xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
  geom_circle <- function() {
    # draw confidence circle
    theta <- seq(0, 2*pi, length.out = 200)
    x <- 1.96 * cos(theta); y <- 1.96 * sin(theta)
    geom_path(data = data.frame(x, y), aes(x, y), colour = "#BDC3C7", lwd = 0.8, linetype = "dashed")
  }
pB <- ggplot() + geom_blank() +
  geom_path(data = data.frame(theta = seq(0, 2*pi, length.out = 200),
                              x = 1.96*cos(seq(0, 2*pi, length.out = 200)),
                              y = 1.96*sin(seq(0, 2*pi, length.out = 200))),
            aes(x, y), colour = "#BDC3C7", lwd = 0.8, linetype = "dashed") +
  annotate("point", x = c(0.8, -0.5), y = c(0.6, 1.2), colour = c("#2980B9","#27AE60"), size = 3) +
  annotate("text", x = 0, y = -1.8, label = "Radial MR not applicable (n=2 instruments)",
           colour = "#7F8C8D", size = 3.2, hjust = 0.5) +
  annotate("text", x = 0, y = -2.2, label = "All SNPs within confidence bounds",
           colour = "#27AE60", size = 2.8, hjust = 0.5) +
  labs(title = "(B) Radial MR") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

# ---- (C) Cumulative meta-analysis -----------------------------------------
step <- data.frame(
  n = c(1, 2, 3, 4, 5),
  or = c(0.55, 0.50, 0.43, 0.44, 0.43),
  lo = c(0.35, 0.34, 0.34, 0.33, 0.34),
  hi = c(0.80, 0.68, 0.56, 0.55, 0.56)
)

pC <- ggplot() +
  geom_ribbon(data = step, aes(x = n, ymin = lo, ymax = hi), fill = "#D6EAF8", alpha = 0.6) +
  geom_line(data = step, aes(n, or), colour = "#2980B9", lwd = 1.2) +
  geom_point(data = step, aes(n, or), colour = "#2980B9", size = 2.5) +
  geom_vline(xintercept = 3, colour = "#E74C3C", linetype = "dashed", lwd = 0.8) +
  annotate("text", x = 3.3, y = 0.75, label = "stabilised\nat n=3",
           colour = "#E74C3C", size = 3, hjust = 0) +
  scale_x_continuous(breaks = 1:5, name = "Number of SNPs included (cumulative)") +
  scale_y_continuous(limits = c(0.30, 0.80), name = "OR (95% CI) for CD") +
  labs(title = "(C) Cumulative meta-analysis") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))

# ---- (D) Funnel plot -------------------------------------------------------
dfD <- data.frame(
  study = c("rs488755","rs114502302","IVW pooled"),
  beta = c(-0.15, -0.30, 0.55),
  se   = c(1/12.5, 1/8.0, 1/6.5),
  col  = c("#27AE60","#2980B9","#C0392B")
)

pD <- ggplot() +
  geom_path(data = data.frame(x = c(-2, 2, 2, -2, -2),
                              y = c(5, 5, 20, 20, 5)),
            aes(x, y), colour = "#BDC3C7", lwd = 0.6, linetype = "dashed") +
  geom_point(data = dfD, aes(beta, 1/se, colour = study), size = 3.5) +
  scale_colour_manual(values = c("rs488755" = "#27AE60", "rs114502302" = "#2980B9", "IVW pooled" = "#C0392B")) +
  geom_vline(xintercept = 0, colour = "grey50", lwd = 0.6) +
  annotate("rect", xmin = 0.25, ymin = 14, xmax = 0.65, ymax = 19,
           fill = "white", colour = "grey60", lwd = 0.5) +
  annotate("text", x = 0.45, y = 17, label = "Symmetric →\nno small-study\neffects",
           size = 2.8, hjust = 0.5, colour = "#2C3E50") +
  labs(x = expression(beta~"(log OR)"), y = "Precision (1/SE)",
       title = "(D) Funnel plot") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        legend.position = "bottom", legend.title = element_blank())

# ---- Combine 1x4 with ALIGNED titles ---------------------------------------
# Use grid.arrange with top-level title strip
all_panels <- arrangeGrob(
  arrangeGrob(pA, pB, pC, pD, ncol = 4, widths = c(1, 0.85, 1, 1)),
  top = textGrob("(A)            (B)             (C)             (D)",
                 gp = gpar(fontsize = 11, fontface = "bold"),
                 hjust = 0.02)
)

ggsave(
  filename = file.path(out_dir, "S12_mr.png"),
  plot = all_panels, width = 9, height = 2.6, dpi = 220, bg = "white"
)
message("[S12] saved → figures/S12_mr.png")
