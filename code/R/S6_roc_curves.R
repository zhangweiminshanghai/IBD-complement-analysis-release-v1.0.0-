# =============================================================================
# S6_roc_curves.R
# Supplementary Figure S6 — ROC curves with Hanley–McNeil 95% CI bands
# *Gut* style: 3 panels A/B/C, single-paragraph caption
# Run: Rscript code/R/S6_roc_curves.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
})

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Hanley–McNeil SE & CI --------------------------------------------------
hanley_se <- function(AUC, n1, n2) {
  # Hanley & McNeil 1982
  Q1 <- AUC / (2 - AUC)
  Q2 <- 2 * AUC^2 / (1 + AUC)
  se <- sqrt((AUC * (1 - AUC) + (n1 - 1) * (Q1 - AUC^2) + (n2 - 1) * (Q2 - AUC^2)) / (n1 * n2))
  se
}

hm_ci <- function(AUC, n1, n2, z = 1.96) {
  se <- hanley_se(AUC, n1, n2)
  data.frame(lo = max(0, AUC - z * se), hi = min(1, AUC + z * se), se = se)
}

# ---- Simulated data (replace with real scores) ------------------------------
set.seed(7)
n_disc_cd <- 73; n_disc_ctrl <- 12
n_val_cd  <- 82; n_val_ctrl <- 15

make_scores <- function(n1, n2, mean_diff) {
  scores <- c(rnorm(n1, mean = mean_diff, sd = 1), rnorm(n2, mean = 0, sd = 1))
  labels <- factor(c(rep("CD", n1), rep("Ctrl", n2)))
  data.frame(scores, labels)
}

d_CFI    <- make_scores(n_disc_cd, n_disc_ctrl, 0.55)
d_PAQR5  <- make_scores(n_disc_cd, n_disc_ctrl, 0.40)
d_KCNE3  <- make_scores(n_disc_cd, n_disc_ctrl, 0.30)
d_3gene  <- make_scores(n_disc_cd, n_disc_ctrl, 0.85)
d_3val   <- make_scores(n_val_cd,  n_val_ctrl,  1.10)

roc_obj <- function(df, ...) {
  roc(df$labels, df$scores, ...)
}

roc_CFI   <- roc_obj(d_CFI,   quiet = TRUE)
roc_PAQR5 <- roc_obj(d_PAQR5, quiet = TRUE)
roc_KCNE3 <- roc_obj(d_KCNE3, quiet = TRUE)
roc_3disc <- roc_obj(d_3gene,  quiet = TRUE)
roc_3val  <- roc_obj(d_3val,   quiet = TRUE)

# ---- AUC + CI table ---------------------------------------------------------
auc_table <- rbind(
  data.frame(gene = "CFI",   AUC = as.numeric(auc(roc_CFI)),   hm_ci(as.numeric(auc(roc_CFI)),   n_disc_cd, n_disc_ctrl)),
  data.frame(gene = "PAQR5", AUC = as.numeric(auc(roc_PAQR5)), hm_ci(as.numeric(auc(roc_PAQR5)), n_disc_cd, n_disc_ctrl)),
  data.frame(gene = "KCNE3", AUC = as.numeric(auc(roc_KCNE3)), hm_ci(as.numeric(auc(roc_KCNE3)), n_disc_cd, n_disc_ctrl)),
  data.frame(gene = "3-gene (disc.)", AUC = as.numeric(auc(roc_3disc)), hm_ci(as.numeric(auc(roc_3disc)), n_disc_cd, n_disc_ctrl)),
  data.frame(gene = "3-gene (val.)",  AUC = as.numeric(auc(roc_3val)),  hm_ci(as.numeric(auc(roc_3val)),  n_val_cd,  n_val_ctrl))
)
print(auc_table)

# Bootstrap CI for validation (per user spec)
boot_auc <- replicate(2000, {
  idx <- sample(seq_along(d_3val$scores), replace = TRUE)
  as.numeric(auc(roc(d_3val$labels[idx], d_3val$scores[idx], quiet = TRUE)))
})
val_boot <- c(AUC = mean(boot_auc),
              lo = quantile(boot_auc, 0.025),
              hi = quantile(boot_auc, 0.975))
message("Bootstrap 3-gene validation AUC: ", sprintf("%.3f [%.3f–%.3f]", val_boot[1], val_boot[2], val_boot[3]))

# ---- Helper: build ROC data.frame with CI band ------------------------------
roc_df <- function(roc_obj) {
  data.frame(FPR = 1 - roc_obj$specificities,
             TPR = roc_obj$sensitivities)
}

# Smooth CI band via Delong variance approximation (use coords)
roc_band <- function(roc_obj, n1, n2, alpha = 0.05) {
  # pointwise approximate SE using Hanley-McNeil binormal assumption
  coords <- coords(roc_obj, x = "all", ret = c("specificity", "sensitivity"))
  fpr <- 1 - coords$specificity
  tpr <- coords$sensitivity
  z <- qnorm(1 - alpha / 2)
  # simplified pointwise SE (Pepe 2003)
  n <- n1 + n2
  se <- sqrt(tpr * (1 - tpr) / n1 + fpr * (1 - fpr) / n2)
  data.frame(FPR = fpr, TPR = tpr, lo = pmax(0, tpr - z * se), hi = pmin(1, tpr + z * se))
}

# ---- Coords -----------------------------------------------------------------
df_CFI   <- roc_band(roc_CFI,   n_disc_cd, n_disc_ctrl)
df_PAQR5 <- roc_band(roc_PAQR5, n_disc_cd, n_disc_ctrl)
df_KCNE3 <- roc_band(roc_KCNE3, n_disc_cd, n_disc_ctrl)
df_3disc <- roc_band(roc_3disc, n_disc_cd, n_disc_ctrl)
df_3val  <- roc_band(roc_3val,  n_val_cd,  n_val_ctrl)

# Step versions for staircase
step_df <- function(df) {
  rbind(data.frame(FPR = 0, TPR = 0, lo = 0, hi = 0),
        data.frame(FPR = df$FPR, TPR = df$TPR, lo = df$lo, hi = df$hi))
}

# ---- Panel A ----------------------------------------------------------------
pA <- ggplot() +
  geom_step(data = step_df(df_CFI),   aes(FPR, TPR), colour = "#E74C3C", lwd = 1.1) +
  geom_ribbon(data = step_df(df_CFI), aes(FPR, ymin = lo, ymax = hi), fill = "#E74C3C", alpha = 0.15) +
  geom_step(data = step_df(df_PAQR5), aes(FPR, TPR), colour = "#2980B9", lwd = 1.1) +
  geom_ribbon(data = step_df(df_PAQR5), aes(FPR, ymin = lo, ymax = hi), fill = "#2980B9", alpha = 0.15) +
  geom_step(data = step_df(df_KCNE3), aes(FPR, TPR), colour = "#27AE60", lwd = 1.1) +
  geom_ribbon(data = step_df(df_KCNE3), aes(FPR, ymin = lo, ymax = hi), fill = "#27AE60", alpha = 0.15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  coord_fixed(ratio = 1) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "(A) Individual genes (GSE16879 discovery)") +
  annotate("text", x = 0.7, y = 0.18, label = sprintf("CFI AUC=%.2f\n[%.2f–%.2f]", auc_table$AUC[1], auc_table$lo[1], auc_table$hi[1]), colour = "#E74C3C", size = 3) +
  annotate("text", x = 0.7, y = 0.05, label = sprintf("PAQR5 AUC=%.2f\n[%.2f–%.2f]", auc_table$AUC[2], auc_table$lo[2], auc_table$hi[2]), colour = "#2980B9", size = 3) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# ---- Panel B ----------------------------------------------------------------
pB <- ggplot() +
  geom_step(data = step_df(df_3disc), aes(FPR, TPR), colour = "#C0392B", lwd = 1.4) +
  geom_ribbon(data = step_df(df_3disc), aes(FPR, ymin = lo, ymax = hi), fill = "#C0392B", alpha = 0.18) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  coord_fixed(ratio = 1) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "(B) 3-gene signature (GSE16879 discovery)") +
  annotate("text", x = 0.55, y = 0.20, label = sprintf("AUC=%.2f\n[%.2f–%.2f]", auc_table$AUC[4], auc_table$lo[4], auc_table$hi[4]), colour = "#C0392B", size = 3.5, fontface = "bold") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# ---- Panel C ----------------------------------------------------------------
pC <- ggplot() +
  geom_step(data = step_df(df_3val), aes(FPR, TPR), colour = "#8E44AD", lwd = 1.4) +
  geom_ribbon(data = step_df(df_3val), aes(FPR, ymin = lo, ymax = hi), fill = "#8E44AD", alpha = 0.18) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  coord_fixed(ratio = 1) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "(C) 3-gene signature (GSE75214 validation)") +
  annotate("text", x = 0.55, y = 0.20, label = sprintf("AUC=%.3f\n[%.3f–%.3f]", val_boot[1], val_boot[2], val_boot[3]), colour = "#8E44AD", size = 3.5, fontface = "bold") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"))

# ---- Combine 1x3 -----------------------------------------------------------
suppressPackageStartupMessages(library(cowplot))
fig <- plot_grid(pA, pB, pC, nrow = 1, rel_widths = c(1, 1, 1),
                 labels = NULL, align = "h")

ggsave(
  filename = file.path(out_dir, "S6_roc.png"),
  plot = fig, width = 9, height = 3.2, dpi = 220, bg = "white"
)
message("[S6] saved → figures/S6_roc.png")
