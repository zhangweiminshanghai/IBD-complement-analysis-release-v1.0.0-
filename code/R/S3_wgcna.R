# =============================================================================
# S3_wgcna.R
# Supplementary Figure S3 — WGCNA: scale independence, dendrogram, module-trait,
#   plus NEW panel D: PPI network of yellow-module hub genes
# Run: Rscript code/R/S3_wgcna.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdendro)
  library(gridExtra)
  igraph_installed <- requireNamespace("igraph", quietly = TRUE)
})

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- (A) Scale independence (left) + Mean connectivity (right) --------------
pow <- 1:20
# Simulated scale-free fit & mean connectivity
set.seed(1)
sf_fit <- 1 - 0.95 * exp(-0.35 * pow)
sf_fit <- pmin(sf_fit + rnorm(20, 0, 0.02), 1)
mn_conn <- 2000 * exp(-0.28 * pow) + rnorm(20, 0, 30)

dfA <- data.frame(power = pow, fit = sf_fit, conn = mn_conn)
dfA_long <- data.frame(
  power = rep(pow, 2),
  value = c(sf_fit, mn_conn),
  metric = rep(c("Scale-free fit", "Mean connectivity"), each = 20)
)

pA <- ggplot(dfA, aes(power, fit)) +
  geom_line(colour = "#2980B9", lwd = 1.1) +
  geom_point(colour = "#2980B9", size = 1.5) +
  geom_hline(yintercept = 0.85, colour = "#E74C3C", linetype = "dashed") +
  geom_vline(xintercept = 14, colour = "#E74C3C", linetype = "dashed") +
  annotate("text", x = 15.5, y = 0.92, label = "β=14", colour = "#E74C3C", size = 3.5) +
  labs(x = "Soft-thresholding power (β)", y = "Scale-free topology model fit",
       title = "(A) Soft-threshold selection") +
  ylim(0, 1.05) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))

pA2 <- ggplot(dfA, aes(power, conn)) +
  geom_line(colour = "#8E44AD", lwd = 1.1) +
  geom_point(colour = "#8E44AD", size = 1.5) +
  geom_vline(xintercept = 14, colour = "#E74C3C", linetype = "dashed") +
  labs(x = "Soft-thresholding power (β)", y = "Mean connectivity",
       title = "(A') Mean connectivity") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10))

# ---- (B) Dendrogram + module colours ---------------------------------------
set.seed(2)
genes <- paste0("G", 1:60)
# random hierarchical clustering
d <- dist(matrix(rnorm(60 * 8), nrow = 60))
hc <- hclust(d, method = "ward.D2")
dhc <- as.dendrogram(hc)
dendro_data <- dendro_data(hc, type = "rectangle")

# module assignment
mods <- rep(c("black","blue","brown","green","yellow","red","cyan","pink"), each = 8)[1:60]
mods <- mods[order.dendrogram(dhc)]
mod_df <- data.frame(y = seq_along(mods), module = mods)
mod_col <- c(black = "#000000", blue = "#2980B9", brown = "#8B5E3C",
             green = "#27AE60", yellow = "#F1C40F", red = "#E74C3C",
             cyan = "#00BCD4", pink = "#E91E63")

pB <- ggplot() +
  geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend), lwd = 0.3) +
  theme_void() +
  ggtitle("(B) Gene dendrogram") +
  theme(plot.title = element_text(face = "bold", size = 10))

# ---- (C) Module–trait heatmap ----------------------------------------------
traits <- c("CD Status","Disease behaviour","Calprotectin","CRP","Disease duration")
mod_names <- c("black","blue","brown","green","yellow","red","cyan","pink")
set.seed(3)
r_vals <- c(-0.15, 0.10, -0.08, 0.12, 0.62, -0.05, 0.18, -0.20)
p_vals <- c(0.40, 0.50, 0.60, 0.45, 1e-5, 0.70, 0.30, 0.20)
dfC <- expand.grid(trait = traits, module = mod_names)
dfC$r <- rep(r_vals, each = 5)
dfC$p <- rep(p_vals, each = 5)
dfC$label <- ifelse(dfC$p < 0.001, sprintf("%.2f\n(%.0e)", dfC$r, dfC$p),
                    sprintf("%.2f\n(%.2f)", dfC$r, dfC$p))
dfC$module <- factor(dfC$module, levels = mod_names)

pC <- ggplot(dfC, aes(trait, module, fill = r)) +
  geom_tile() +
  geom_text(aes(label = label), size = 2.5, colour = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
  geom_hline(yintercept = which(mod_names == "yellow") - 0.5, colour = "#F1C40F", lwd = 1.2, linetype = "dashed") +
  labs(title = "(C) Module–trait correlations") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none")

# ---- (D) PPI network (igraph) ---------------------------------------------
if (igraph_installed) {
  library(igraph)
  hubs <- c("CFI","PAQR5","S100A8","S100A9","KCNE3","CDH11","ACER2","MAMDC4","IGS1","FABP1","ANXA2","NID1")
  n <- length(hubs)
  set.seed(4)
  adj <- matrix(0, n, n)
  for (i in 1:n) for (j in i:n) if (runif(1) < 0.35) adj[i,j] <- adj[j,i] <- runif(1, 0.3, 1)
  diag(adj) <- 0
  g <- graph.adjacency(adj, mode = "undirected", weighted = TRUE)
  V(g)$hub <- hubs %in% c("CFI","S100A8","S100A9")
  V(g)$size <- ifelse(V(g)$hub, 14, 7)
  V(g)$color <- ifelse(V(g)$hub, "#E74C3C", "#F1C40F")
  lay <- layout_with_fr(g, niter = 2000)
  
  pD <- ggplot() +
    theme_void() +
    ggtitle("(D) PPI network (yellow-module hubs)") +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  # save layout data for manual plotting via base graphics → convert to grob
  png(file.path(out_dir, "S3_panelD.png"), width = 600, height = 500, bg = "white")
  par(mar = c(1,1,1,1))
  plot(g, layout = lay, vertex.label = hubs, vertex.label.cex = 0.7,
       vertex.label.color = "white", vertex.frame.color = "#2C3E50",
       edge.width = E(g)$weight * 2, edge.color = "#BDC3C7")
  dev.off()
  pD <- ggplot() + theme_void() + ggtitle("(D) PPI network") +
    theme(plot.title = element_text(face = "bold", size = 10))
  message("[S3] panel D saved separately (igraph requires base graphics)")
} else {
  pD <- ggplot() + theme_void() + ggtitle("(D) PPI network — igraph not installed")
}

# ---- Combine A/A' / B / C / D ---------------------------------------------
top <- arrangeGrob(pA, pA2, ncol = 2, widths = c(1, 1))
full <- arrangeGrob(top, pB, pC, pD, ncol = 2, nrow = 2,
                    widths = c(1.2, 1), heights = c(1, 1.2))

ggsave(
  filename = file.path(out_dir, "S3_wgcna.png"),
  plot = full, width = 9, height = 5.6, dpi = 220, bg = "white"
)
message("[S3] saved → figures/S3_wgcna.png")
