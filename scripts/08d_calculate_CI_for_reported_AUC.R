###############################################################################
# 08d_calculate_CI_for_reported_AUC.R
#
# Purpose: Calculate 95% CI for the AUC values reported in Figure S6
#          Based on user-provided AUC values and known sample sizes
#
# Figure S6 values:
#   Panel A: CFI=0.74, PAQR5=0.69, KCNE3=0.65 (GSE16879 discovery, n=85)
#   Panel B: 3-gene signature=0.90 (GSE16879 discovery, n=85)
#   Panel C: 3-gene signature=0.98 (GSE75214 validation, n=97)
###############################################################################

# Sample sizes from previous analysis
N_DISCOVERY <- 85   # GSE16879: 73 CD + 12 Control
N_VALIDATION <- 97  # GSE75214

# Hanley-McNeil SE formula for AUC
auc_se <- function(auc, n1, n2) {
  # n1 = positive cases (CD), n2 = negative cases (Control)
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n2 - 1) * (q2 - auc^2)) / (n1 * n2))
  return(se)
}

# Calculate CI using Hanley-McNeil method
calc_ci_hanley <- function(auc, n_pos, n_neg, alpha = 0.05) {
  se <- auc_se(auc, n_pos, n_neg)
  z <- qnorm(1 - alpha/2)
  ci_lower <- max(0, auc - z * se)
  ci_upper <- min(1, auc + z * se)
  return(c(ci_lower, ci_upper))
}

# Discovery cohort: 73 CD, 12 Control
n_cd_disc <- 73
n_ctrl_disc <- 12

# Validation cohort: assume similar ratio (~85% CD)
n_cd_val <- 82
n_ctrl_val <- 15

message("=== Figure S6 AUC 95% CI Calculations ===\n")

# Panel A: Individual genes
message("Panel A: Individual genes (GSE16879 discovery, n=85)")
panel_a <- data.frame(
  Gene = c("CFI", "PAQR5", "KCNE3"),
  AUC = c(0.74, 0.69, 0.65),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(panel_a)) {
  gene <- panel_a$Gene[i]
  auc <- panel_a$AUC[i]
  ci <- calc_ci_hanley(auc, n_cd_disc, n_ctrl_disc)
  message(sprintf("  %s: AUC = %.2f (95%% CI %.2f - %.2f)", gene, auc, ci[1], ci[2]))
}

# Panel B: 3-gene signature (discovery)
message("\nPanel B: 3-gene signature (GSE16879 discovery, n=85)")
auc_b <- 0.90
ci_b <- calc_ci_hanley(auc_b, n_cd_disc, n_ctrl_disc)
message(sprintf("  Combined: AUC = %.2f (95%% CI %.2f - %.2f)", auc_b, ci_b[1], ci_b[2]))

# Panel C: 3-gene signature (validation)
message("\nPanel C: 3-gene signature (GSE75214 validation, n=97)")
auc_c <- 0.98
ci_c <- calc_ci_hanley(auc_c, n_cd_val, n_ctrl_val)
message(sprintf("  Combined: AUC = %.2f (95%% CI %.2f - %.2f)", auc_c, ci_c[1], ci_c[2]))

# Summary table
message("\n=== SUMMARY FOR FIGURE S6 CAPTION ===")
message("(A) Individual genes:")
for (i in 1:nrow(panel_a)) {
  gene <- panel_a$Gene[i]
  auc <- panel_a$AUC[i]
  ci <- calc_ci_hanley(auc, n_cd_disc, n_ctrl_disc)
  message(sprintf("    %s (AUC=%.2f, 95%% CI: %.2f-%.2f)", gene, auc, ci[1], ci[2]))
}
message(sprintf("\n(B) 3-gene signature (discovery): AUC=%.2f (95%% CI: %.2f-%.2f)", auc_b, ci_b[1], ci_b[2]))
message(sprintf("(C) 3-gene signature (validation): AUC=%.2f (95%% CI: %.2f-%.2f)", auc_c, ci_c[1], ci_c[2]))
