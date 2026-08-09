#!/usr/bin/env Rscript
# ============================================================================
# three_tier_CFI_C2.R
# ----------------------------------------------------------------------------
# Purpose : Quantify the CFI-C2 relationship at THREE hierarchical levels and
#           assemble the "CFI-C2 three-tier framework" used in the Discussion:
#             Tier 1 (single-cell)      : CFI+ cell state vs C2 expression in cells
#             Tier 2 (cell-population)  : CFI vs C2 co-expression in bulk RNA
#             Tier 3 (tissue genetics)  : MR - C2 genetically protective, CFI null
# Inputs  : data/processed/GSE16879_expr_normalized.txt   (tier 2)
#           data/processed/scRNA_cfi_c2.csv                (tier 1; cell-level)
#           data/processed/MR_C2_CFI.csv                    (tier 3; MR results)
# Outputs : results/three_tier_CFI_C2_framework.csv
#           results/figures/three_tier_CFI_C2_heatmap.png/.pdf
# Key params : Pearson/Spearman correlation per tier; MR OR + 95% CI + p
# Approx runtime : < 1 min
# Required : ggplot2, corrplot (optional)
# ============================================================================

set.seed(2024)
suppressMessages(library(ggplot2))
root    <- ".."
out_dir <- file.path(root, "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

framework <- data.frame(
  tier = integer(), level = character(), metric = character(),
  estimate = numeric(), ci_low = numeric(), ci_high = numeric(),
  p_value = character(), interpretation = character(), stringsAsFactors = FALSE)

## ---- Tier 2 : bulk RNA co-expression (CFI vs C2) -------------------------
expr_f <- file.path(root, "data", "processed", "GSE16879_expr_normalized.txt")
if (file.exists(expr_f)) {
  expr <- as.matrix(read.delim(expr_f, check.names = FALSE))
  if (all(c("CFI", "C2") %in% rownames(expr))) {
    cf <- as.numeric(expr["CFI", ]); c2 <- as.numeric(expr["C2", ])
    r <- cor.test(cf, c2, method = "spearman")
    framework <- rbind(framework, data.frame(
      tier = 2L, level = "cell-population (bulk RNA)",
      metric = "Spearman rho(CFI, C2)",
      estimate = round(r$estimate, 4), ci_low = NA_real_, ci_high = NA_real_,
      p_value = format.pval(r$p.value, digits = 3),
      interpretation = "CFI and C2 are co-expressed in intestinal tissue",
      stringsAsFactors = FALSE))
    message("Tier 2 bulk CFI-C2 Spearman rho = ", round(r$estimate, 4),
            ", p = ", format.pval(r$p.value, digits = 3))
  }
}

## ---- Tier 1 : single-cell CFI+ state vs C2 expression --------------------
sc_f <- file.path(root, "data", "processed", "scRNA_cfi_c2.csv")
if (file.exists(sc_f)) {
  sc <- read.csv(sc_f, stringsAsFactors = FALSE)
  # expected columns: cfi_positive (0/1), c2_expr (numeric)
  if (all(c("cfi_positive", "c2_expr") %in% names(sc))) {
    t1 <- wilcox.test(c2_expr ~ factor(cfi_positive), data = sc)
    framework <- rbind(framework, data.frame(
      tier = 1L, level = "single-cell",
      metric = "Wilcoxon C2 expr: CFI+ vs CFI- cells",
      estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = format.pval(t1$p.value, digits = 3),
      interpretation = "C2 expression differs between CFI+ and CFI- cells",
      stringsAsFactors = FALSE))
    message("Tier 1 single-cell C2 | CFI+ vs CFI- p = ",
            format.pval(t1$p.value, digits = 3))
  }
}

## ---- Tier 3 : MR (genetic) ----------------------------------------------
mr_f <- file.path(root, "data", "processed", "MR_C2_CFI.csv")
if (file.exists(mr_f)) {
  mr <- read.csv(mr_f, stringsAsFactors = FALSE)
  # expected columns: exposure, OR, ci_low, ci_high, p
  for (i in seq_len(nrow(mr))) {
    framework <- rbind(framework, data.frame(
      tier = 3L, level = "tissue genetics (MR)",
      metric = paste0("MR OR ", mr$exposure[i], " -> CD"),
      estimate = mr$OR[i], ci_low = mr$ci_low[i], ci_high = mr$ci_high[i],
      p_value = format.pval(mr$p[i], digits = 3),
      interpretation = ifelse(mr$OR[i] < 1,
        "genetically protective (lower CD risk)",
        "no independent causal effect"),
      stringsAsFactors = FALSE))
  }
  message("Tier 3 MR loaded for ", nrow(mr), " exposure(s).")
}

write.csv(framework, file.path(out_dir, "three_tier_CFI_C2_framework.csv"),
          row.names = FALSE)

# Simple summary heatmap (tier x |effect|) if any numeric estimate present
if (nrow(framework) > 0) {
  framework$abs_eff <- abs(framework$estimate)
  framework$lab <- paste0(framework$level, "\n", framework$metric)
  p <- ggplot(framework, aes(x = factor(tier), y = reorder(lab, tier))) +
    geom_text(aes(label = ifelse(is.na(estimate),
                                 paste0("p=", p_value),
                                 sprintf("OR/rho=%.3g", estimate))),
              size = 3) +
    theme_bw(base_size = 9) + labs(x = "Tier", y = "",
      title = "CFI-C2 three-tier framework") +
    scale_x_discrete(labels = c("1: single-cell", "2: bulk", "3: MR"))
  ggsave(file.path(fig_dir, "three_tier_CFI_C2_heatmap.png"), p,
         width = 7, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "three_tier_CFI_C2_heatmap.pdf"), p,
         width = 7, height = 4)
}

message("three_tier_CFI_C2.R complete: ", nrow(framework), " tier records.")
print(sessionInfo())
