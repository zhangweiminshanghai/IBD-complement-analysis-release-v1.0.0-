###############################################################################
# 06_MR_TwoSampleMR.R
#
# Purpose      : Two-sample Mendelian randomisation testing whether genetically
#                predicted complement protein levels are causally related to
#                Crohn's disease (CD):
#                  * C2  -> CD : protective causal effect (OR = 0.43 per SD,
#                                p = 2.2e-6, inverse-variance weighted)
#                  * CFI -> CD : no independent causal effect (p = 0.73)
#                  * CD  -> CFI (reverse direction): compensatory up-regulation
#                Sensitivity analyses: MR-Egger intercept (pleiotropy), Cochran's Q
#                (heterogeneity), weighted median, leave-one-out.
#
# Inputs       : ONLINE  (default) - OpenGWAS via TwoSampleMR/ieugwasr.
#                    Requires an OpenGWAS JWT in the environment variable
#                    OPENGWAS_JWT (see ieugwasr::get_opengwas_jwt()).
#                OFFLINE (fallback) - pre-extracted harmonised instrument tables:
#                    data/processed/mr/C2_exposure_instruments.csv
#                    data/processed/mr/CFI_exposure_instruments.csv
#                    data/processed/mr/CD_outcome_for_C2.csv
#                    data/processed/mr/CD_outcome_for_CFI.csv
#                    data/processed/mr/CFI_outcome_for_CD.csv   (reverse MR, optional)
#                  Required columns (TwoSampleMR naming):
#                    SNP, beta.exposure/outcome, se.exposure/outcome,
#                    effect_allele.*, other_allele.*, eaf.*, pval.*
#
# Outputs      : results/06_MR/mr_results.csv              (all methods, both directions)
#                results/06_MR/mr_odds_ratios.csv
#                results/06_MR/mr_heterogeneity.csv
#                results/06_MR/mr_pleiotropy.csv
#                results/06_MR/mr_single_snp.csv
#                results/06_MR/mr_leaveoneout.csv
#                results/06_MR/harmonised_<exposure>_<outcome>.csv
#                results/figures/MR_forest_plot.{png,pdf}
#                results/figures/MR_scatter_<exposure>_CD.{png,pdf}
#                results/figures/MR_leaveoneout_<exposure>_CD.{png,pdf}
#
# Figure/Table : Figure 6A-C ; Table 3 (MR estimates and sensitivity analyses)
#
# Key params   : Instrument selection : p < 5e-8 (relaxed to 5e-6 for pQTL when
#                    fewer than MIN_SNPS instruments survive), clump r2 = 0.001,
#                    window = 10,000 kb, EUR reference panel
#                Primary method       : inverse-variance weighted (multiplicative
#                    random effects when nsnp >= 3, Wald ratio when nsnp == 1)
#                Effect scale         : OR of CD per 1-SD higher protein level
#                set.seed(2024)
#
# Runtime      : ~5-15 min online (API dependent); < 1 min offline
#
# Packages     : TwoSampleMR (>= 0.5.11), ieugwasr (>= 1.0), ggplot2 (>= 3.4),
#                data.table (>= 1.14); optional MendelianRandomization (>= 0.9)
#
# Author       : IBD complement project
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

set.seed(2024)

## ---------------------------------------------------------------------------
## 0. Configuration
## ---------------------------------------------------------------------------
find_repo_root <- function() {
  cand <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:4) {
    if (dir.exists(file.path(cand, "data")) && dir.exists(file.path(cand, "scripts"))) return(cand)
    cand <- normalizePath(file.path(cand, ".."), winslash = "/", mustWork = FALSE)
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
ROOT     <- find_repo_root()
MR_DIR   <- file.path(ROOT, "data", "processed", "mr")
RES_DIR  <- file.path(ROOT, "results", "06_MR")
FIG_DIR  <- file.path(ROOT, "results", "figures")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# OpenGWAS study identifiers -------------------------------------------------
# Exposures: circulating complement proteins (pQTL, Sun et al. 2018 INTERVAL)
#   C2  : prot-a-435   (complement C2)        - fallback eQTL eqtl-a-ENSG00000166278
#   CFI : prot-a-523   (complement factor I)  - fallback eQTL eqtl-a-ENSG00000205403
# Outcome: Crohn's disease
#   ieu-a-30           (de Lange et al. 2017, 12,194 cases / 28,072 controls)
#   finn-b-K11_CROHN   (FinnGen R5) used as a replication outcome
EXPOSURES <- list(
  C2  = list(id = "prot-a-435", fallback = "eqtl-a-ENSG00000166278", label = "C2 (per SD)"),
  CFI = list(id = "prot-a-523", fallback = "eqtl-a-ENSG00000205403", label = "CFI (per SD)")
)
OUTCOME_CD      <- "ieu-a-30"
OUTCOME_CD_REPL <- "finn-b-K11_CROHN"
CFI_OUTCOME_ID  <- "prot-a-523"    # for the reverse direction CD -> CFI

P_THRESHOLD     <- 5e-8
P_RELAXED       <- 5e-6
CLUMP_R2        <- 0.001
CLUMP_KB        <- 10000
MIN_SNPS        <- 2

have_tsmr <- requireNamespace("TwoSampleMR", quietly = TRUE)
if (have_tsmr) suppressPackageStartupMessages(library(TwoSampleMR))
online <- have_tsmr && nzchar(Sys.getenv("OPENGWAS_JWT"))
message("TwoSampleMR available: ", have_tsmr, " | OpenGWAS token set: ",
        nzchar(Sys.getenv("OPENGWAS_JWT")))
if (!online) {
  message("-> running in OFFLINE mode: reading pre-extracted instruments from ", MR_DIR)
}

## ---------------------------------------------------------------------------
## 1. Manual IVW / Wald-ratio implementation (used offline and as a cross-check)
## ---------------------------------------------------------------------------
#' Fixed/multiplicative random-effects inverse-variance weighted MR.
#' @param bx,bxse SNP-exposure effects and standard errors
#' @param by,byse SNP-outcome effects and standard errors
ivw_manual <- function(bx, bxse, by, byse) {
  keep <- stats::complete.cases(bx, bxse, by, byse) & byse > 0
  bx <- bx[keep]; bxse <- bxse[keep]; by <- by[keep]; byse <- byse[keep]
  n  <- length(bx)
  if (n == 0) return(NULL)
  if (n == 1) {
    b  <- by / bx
    se <- abs(byse / bx)                                   # first-order Wald ratio SE
    return(data.frame(method = "Wald ratio", nsnp = 1, b = b, se = se,
                      pval = 2 * stats::pnorm(-abs(b / se))))
  }
  w   <- 1 / byse^2
  b   <- sum(w * bx * by) / sum(w * bx^2)
  se  <- sqrt(1 / sum(w * bx^2))
  # Cochran's Q and multiplicative random-effects inflation
  Q   <- sum(w * (by - b * bx)^2)
  Qdf <- n - 1
  phi <- max(1, Q / Qdf)
  se_re <- se * sqrt(phi)
  data.frame(method = "Inverse variance weighted", nsnp = n, b = b, se = se_re,
             pval = 2 * stats::pnorm(-abs(b / se_re)),
             Q = Q, Q_df = Qdf, Q_pval = stats::pchisq(Q, Qdf, lower.tail = FALSE))
}

#' Simple MR-Egger for offline mode (intercept = directional pleiotropy).
egger_manual <- function(bx, bxse, by, byse) {
  keep <- stats::complete.cases(bx, bxse, by, byse) & byse > 0
  bx <- bx[keep]; by <- by[keep]; byse <- byse[keep]
  if (length(bx) < 3) return(NULL)
  # orient exposure effects positive (standard MR-Egger convention)
  sgn <- sign(bx); sgn[sgn == 0] <- 1
  fit <- stats::lm(I(by * sgn) ~ I(bx * sgn), weights = 1 / byse^2)
  s   <- summary(fit)$coefficients
  data.frame(method = "MR Egger", nsnp = length(bx),
             b = s[2, 1], se = s[2, 2], pval = s[2, 4],
             egger_intercept = s[1, 1], intercept_se = s[1, 2], intercept_pval = s[1, 4])
}

add_or <- function(df) {
  df$OR       <- exp(df$b)
  df$OR_lci95 <- exp(df$b - 1.96 * df$se)
  df$OR_uci95 <- exp(df$b + 1.96 * df$se)
  df
}

## ---------------------------------------------------------------------------
## 2. Instrument extraction (online) or file loading (offline)
## ---------------------------------------------------------------------------
get_exposure_online <- function(cfg, name) {
  inst <- try(TwoSampleMR::extract_instruments(outcomes = cfg$id, p1 = P_THRESHOLD,
                                               clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB),
              silent = TRUE)
  if (inherits(inst, "try-error") || is.null(inst) || nrow(inst) < MIN_SNPS) {
    message("  relaxing threshold to ", P_RELAXED, " for ", name)
    inst <- try(TwoSampleMR::extract_instruments(outcomes = cfg$id, p1 = P_RELAXED,
                                                 clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB),
                silent = TRUE)
  }
  if (inherits(inst, "try-error") || is.null(inst) || nrow(inst) < 1) {
    message("  falling back to the eQTL dataset ", cfg$fallback)
    inst <- try(TwoSampleMR::extract_instruments(outcomes = cfg$fallback, p1 = P_RELAXED,
                                                 clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB),
                silent = TRUE)
  }
  if (inherits(inst, "try-error") || is.null(inst)) return(NULL)
  inst$exposure <- name
  inst
}

read_offline <- function(file) {
  p <- file.path(MR_DIR, file)
  if (!file.exists(p)) return(NULL)
  data.table::fread(p, data.table = FALSE)
}

## ---------------------------------------------------------------------------
## 3. Run MR for each exposure -> CD
## ---------------------------------------------------------------------------
all_res  <- list()
all_het  <- list()
all_plt  <- list()
all_snp  <- list()
all_loo  <- list()

for (nm in names(EXPOSURES)) {
  message("\n=== MR: ", nm, " -> Crohn's disease ===")
  cfg <- EXPOSURES[[nm]]
  dat <- NULL

  if (online) {
    exp_dat <- get_exposure_online(cfg, nm)
    if (!is.null(exp_dat) && nrow(exp_dat)) {
      message("  instruments: ", nrow(exp_dat), " SNPs")
      out_dat <- try(TwoSampleMR::extract_outcome_data(snps = exp_dat$SNP,
                                                       outcomes = OUTCOME_CD,
                                                       proxies = TRUE), silent = TRUE)
      if (inherits(out_dat, "try-error") || is.null(out_dat)) {
        out_dat <- try(TwoSampleMR::extract_outcome_data(snps = exp_dat$SNP,
                                                         outcomes = OUTCOME_CD_REPL,
                                                         proxies = TRUE), silent = TRUE)
      }
      if (!inherits(out_dat, "try-error") && !is.null(out_dat) && nrow(out_dat)) {
        dat <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = 2)
        dat <- dat[dat$mr_keep, , drop = FALSE]
      }
    }
  }

  if (is.null(dat) || !nrow(dat)) {
    e <- read_offline(paste0(nm, "_exposure_instruments.csv"))
    o <- read_offline(paste0("CD_outcome_for_", nm, ".csv"))
    if (is.null(e) || is.null(o)) {
      warning("No data for ", nm, " (online failed and offline files missing) - skipped.")
      next
    }
    if (have_tsmr && all(c("effect_allele.exposure") %in% colnames(e))) {
      dat <- TwoSampleMR::harmonise_data(e, o, action = 2)
      dat <- dat[dat$mr_keep, , drop = FALSE]
    } else {
      dat <- merge(e, o, by = "SNP")
    }
    message("  offline instruments: ", nrow(dat), " SNPs")
  }

  data.table::fwrite(dat, file.path(RES_DIR, paste0("harmonised_", nm, "_CD.csv")))

  if (have_tsmr && all(c("id.exposure", "id.outcome") %in% colnames(dat))) {
    res <- TwoSampleMR::mr(dat, method_list = c(
      "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode",
      "mr_wald_ratio"))
    res$exposure <- nm; res$outcome <- "Crohn's disease"
    all_res[[nm]] <- add_or(res[, c("exposure", "outcome", "method", "nsnp", "b", "se", "pval")])

    het <- try(TwoSampleMR::mr_heterogeneity(dat), silent = TRUE)
    if (!inherits(het, "try-error") && !is.null(het)) { het$exposure <- nm; all_het[[nm]] <- het }
    plt <- try(TwoSampleMR::mr_pleiotropy_test(dat), silent = TRUE)
    if (!inherits(plt, "try-error") && !is.null(plt)) { plt$exposure <- nm; all_plt[[nm]] <- plt }
    ss  <- try(TwoSampleMR::mr_singlesnp(dat), silent = TRUE)
    if (!inherits(ss, "try-error") && !is.null(ss)) { ss$exposure <- nm; all_snp[[nm]] <- ss }
    loo <- try(TwoSampleMR::mr_leaveoneout(dat), silent = TRUE)
    if (!inherits(loo, "try-error") && !is.null(loo)) { loo$exposure <- nm; all_loo[[nm]] <- loo }

    # figures
    ps <- try(TwoSampleMR::mr_scatter_plot(res, dat), silent = TRUE)
    if (!inherits(ps, "try-error")) {
      ggsave(file.path(FIG_DIR, paste0("MR_scatter_", nm, "_CD.png")), ps[[1]],
             width = 7, height = 6, dpi = 300)
      ggsave(file.path(FIG_DIR, paste0("MR_scatter_", nm, "_CD.pdf")), ps[[1]],
             width = 7, height = 6)
    }
    if (!is.null(all_loo[[nm]])) {
      pl <- try(TwoSampleMR::mr_leaveoneout_plot(all_loo[[nm]]), silent = TRUE)
      if (!inherits(pl, "try-error")) {
        ggsave(file.path(FIG_DIR, paste0("MR_leaveoneout_", nm, "_CD.png")), pl[[1]],
               width = 7, height = max(4, 0.25 * nrow(all_loo[[nm]])), dpi = 300, limitsize = FALSE)
        ggsave(file.path(FIG_DIR, paste0("MR_leaveoneout_", nm, "_CD.pdf")), pl[[1]],
               width = 7, height = max(4, 0.25 * nrow(all_loo[[nm]])), limitsize = FALSE)
      }
    }
  } else {
    # offline / manual computation
    bx <- dat[["beta.exposure"]]; bxse <- dat[["se.exposure"]]
    by <- dat[["beta.outcome"]];  byse <- dat[["se.outcome"]]
    r  <- ivw_manual(bx, bxse, by, byse)
    e2 <- egger_manual(bx, bxse, by, byse)
    res <- do.call(rbind, lapply(list(r, e2), function(z) if (is.null(z)) NULL else
      data.frame(exposure = nm, outcome = "Crohn's disease",
                 method = z$method, nsnp = z$nsnp, b = z$b, se = z$se, pval = z$pval)))
    all_res[[nm]] <- add_or(res)
    if (!is.null(r) && "Q" %in% names(r)) {
      all_het[[nm]] <- data.frame(exposure = nm, method = "IVW", Q = r$Q,
                                  Q_df = r$Q_df, Q_pval = r$Q_pval)
    }
    if (!is.null(e2)) {
      all_plt[[nm]] <- data.frame(exposure = nm, egger_intercept = e2$egger_intercept,
                                  se = e2$intercept_se, pval = e2$intercept_pval)
    }
  }

  r_ivw <- all_res[[nm]][all_res[[nm]]$method %in%
                           c("Inverse variance weighted", "Wald ratio"), ]
  if (nrow(r_ivw)) {
    message(sprintf("  %s -> CD: OR = %.3f (95%% CI %.3f-%.3f), P = %.3g, nSNP = %d",
                    nm, r_ivw$OR[1], r_ivw$OR_lci95[1], r_ivw$OR_uci95[1],
                    r_ivw$pval[1], r_ivw$nsnp[1]))
  }
}

## ---------------------------------------------------------------------------
## 4. Reverse direction: CD -> CFI (compensatory up-regulation)
## ---------------------------------------------------------------------------
message("\n=== Reverse MR: Crohn's disease -> CFI ===")
rev_res <- NULL
if (online) {
  exp_cd <- try(TwoSampleMR::extract_instruments(outcomes = OUTCOME_CD, p1 = P_THRESHOLD,
                                                 clump = TRUE, r2 = CLUMP_R2, kb = CLUMP_KB),
                silent = TRUE)
  if (!inherits(exp_cd, "try-error") && !is.null(exp_cd) && nrow(exp_cd)) {
    out_cfi <- try(TwoSampleMR::extract_outcome_data(snps = exp_cd$SNP,
                                                     outcomes = CFI_OUTCOME_ID), silent = TRUE)
    if (!inherits(out_cfi, "try-error") && !is.null(out_cfi) && nrow(out_cfi)) {
      dat_rev <- TwoSampleMR::harmonise_data(exp_cd, out_cfi, action = 2)
      dat_rev <- dat_rev[dat_rev$mr_keep, , drop = FALSE]
      data.table::fwrite(dat_rev, file.path(RES_DIR, "harmonised_CD_CFI.csv"))
      rr <- TwoSampleMR::mr(dat_rev, method_list = c("mr_ivw", "mr_weighted_median"))
      rev_res <- add_or(data.frame(exposure = "Crohn's disease", outcome = "CFI",
                                   method = rr$method, nsnp = rr$nsnp,
                                   b = rr$b, se = rr$se, pval = rr$pval))
    }
  }
}
if (is.null(rev_res)) {
  e <- read_offline("CD_exposure_instruments.csv")
  o <- read_offline("CFI_outcome_for_CD.csv")
  if (!is.null(e) && !is.null(o)) {
    d <- merge(e, o, by = "SNP")
    r <- ivw_manual(d$beta.exposure, d$se.exposure, d$beta.outcome, d$se.outcome)
    if (!is.null(r)) {
      rev_res <- add_or(data.frame(exposure = "Crohn's disease", outcome = "CFI",
                                   method = r$method, nsnp = r$nsnp, b = r$b,
                                   se = r$se, pval = r$pval))
    }
  } else {
    message("  reverse-direction data not available - skipped.")
  }
}
if (!is.null(rev_res)) {
  all_res[["CD_to_CFI"]] <- rev_res
  message(sprintf("  CD -> CFI: beta = %.4f, P = %.3g", rev_res$b[1], rev_res$pval[1]))
}

## ---------------------------------------------------------------------------
## 5. Write result tables
## ---------------------------------------------------------------------------
if (!length(all_res)) stop("No MR results were produced - check inputs / OpenGWAS access.")
res_all <- do.call(rbind, all_res)
rownames(res_all) <- NULL
data.table::fwrite(res_all, file.path(RES_DIR, "mr_results.csv"))
data.table::fwrite(res_all[res_all$method %in% c("Inverse variance weighted", "Wald ratio"), ],
                   file.path(RES_DIR, "mr_odds_ratios.csv"))
if (length(all_het)) data.table::fwrite(do.call(rbind, all_het),
                                        file.path(RES_DIR, "mr_heterogeneity.csv"))
if (length(all_plt)) data.table::fwrite(do.call(rbind, all_plt),
                                        file.path(RES_DIR, "mr_pleiotropy.csv"))
if (length(all_snp)) data.table::fwrite(do.call(rbind, all_snp),
                                        file.path(RES_DIR, "mr_single_snp.csv"))
if (length(all_loo)) data.table::fwrite(do.call(rbind, all_loo),
                                        file.path(RES_DIR, "mr_leaveoneout.csv"))
print(res_all)

## ---------------------------------------------------------------------------
## 6. Forest plot of the primary estimates
## ---------------------------------------------------------------------------
fp <- res_all[res_all$method %in% c("Inverse variance weighted", "Wald ratio") &
                res_all$outcome == "Crohn's disease", , drop = FALSE]
if (nrow(fp)) {
  fp$Label <- sprintf("%s -> CD (%d SNPs)", fp$exposure, fp$nsnp)
  fp$Annot <- sprintf("OR = %.2f (%.2f-%.2f), P = %.2g",
                      fp$OR, fp$OR_lci95, fp$OR_uci95, fp$pval)
  fp$Signif <- ifelse(fp$pval < 0.05, "P < 0.05", "n.s.")

  p_forest <- ggplot(fp, aes(x = OR, y = Label, colour = Signif)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = OR_lci95, xmax = OR_uci95), orientation = "y",
                  width = 0.12, linewidth = 0.7) +
    geom_point(size = 3.2) +
    geom_text(aes(label = Annot), vjust = -1.1, size = 3.2, show.legend = FALSE) +
    scale_x_continuous(trans = "log10") +
    scale_colour_manual(values = c("P < 0.05" = "#2E86AB", "n.s." = "grey45")) +
    labs(title = "Mendelian randomisation: complement proteins and Crohn's disease",
         subtitle = "Inverse-variance weighted estimates, OR per 1-SD higher protein level",
         x = "Odds ratio (log scale)", y = NULL, colour = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  ggsave(file.path(FIG_DIR, "MR_forest_plot.png"), p_forest, width = 9, height = 4.5, dpi = 300)
  ggsave(file.path(FIG_DIR, "MR_forest_plot.pdf"), p_forest, width = 9, height = 4.5)
}

message("\nReference values reported in the manuscript:")
message("  C2  -> CD : OR = 0.43 per SD, P = 2.2e-6 (protective, IVW)")
message("  CFI -> CD : P = 0.73 (no independent causal effect; compensatory marker)")
message("\nDone. Tables -> ", RES_DIR, " ; figures -> ", FIG_DIR)

## ---------------------------------------------------------------------------
## 7. Session information
## ---------------------------------------------------------------------------
sessionInfo()
