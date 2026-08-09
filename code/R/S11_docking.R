# =============================================================================
# S11_docking.R
# Supplementary Figure S11 — Nafamostat binding interactions (C-panel fixed)
# *Gut* style: strong binders grouped, weak/non-specific grouped, no overlap
# Run: Rscript code/R/S11_docking.R
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Binding free energy data (kcal/mol) -----------------------------------
df <- data.frame(
  target = factor(c("TNFSF10","SGS","TYRO3","TIRAP","TNKS4","TRAF2"),
                  levels = c("TNFSF10","SGS","TYRO3","TIRAP","TNKS4","TRAF2")),
  dG      = c(1.2, 1.5, 8.9, 9.7, 9.9, 9.6),
  group   = c("weak","weak","strong","strong","strong","strong")
)

# reorder: weak first (bottom), strong on top
df$target <- factor(df$target, levels = c("TNFSF10","SGS","TYRO3","TRAF2","TIRAP","TNKS4"))

# colour mapping
col_map <- c(weak = "#BDC3C7", strong = c(TYRO3 = "#8E44AD", TRAF2 = "#2980B9",
                                            TIRAP = "#27AE60", TNKS4 = "#C0392B"))

p <- ggplot(df, aes(x = dG, y = target, fill = group)) +
  geom_bar(stat = "identity", width = 0.55, colour = "#2C3E50", size = 0.3) +
  scale_fill_manual(values = c("weak" = "#BDC3C7",
                                "strong_TYRO3" = "#8E44AD",
                                "strong_TRAF2" = "#2980B9",
                                "strong_TIRAP" = "#27AE60",
                                "strong_TNKS4" = "#C0392B")) +
  geom_text(aes(label = sprintf("%.1f", dG)), hjust = -0.3, size = 3.5, colour = "#2C3E50") +
  geom_vline(xintercept = 7, colour = "#E74C3C", linetype = "dashed", lwd = 0.8) +
  annotate("text", x = 7.5, y = 6.3, label = "strong binders →", colour = "#E74C3C",
           size = 3.2, hjust = 0) +
  # annotation for C1r/C1s related
  annotate("text", x = 9.8, y = 3.5, label = "C1r/C1s related", colour = "#2980B9",
           size = 3, fontface = "italic") +
  annotate("text", x = 9.8, y = 2.5, label = "Weak/non-specific", colour = "#7F8C8D",
           size = 3, fontface = "italic") +
  scale_x_continuous(limits = c(0, 12), expand = c(0.02, 0.02),
                     name = expression(paste("Binding free energy  ", Delta, "G (kcal/mol)"))) +
  labs(y = "", title = "(C) Nafamostat binding interactions") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "none",
        axis.text.y = element_text(size = 9, colour = "#2C3E50"),
        panel.grid.minor = element_blank())

ggsave(
  filename = file.path(out_dir, "S11_docking.png"),
  plot = p, width = 9, height = 3.5, dpi = 220, bg = "white"
)
message("[S11] saved → figures/S11_docking.png")
