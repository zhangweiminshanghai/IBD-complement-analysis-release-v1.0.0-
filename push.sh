#!/usr/bin/env bash
# =============================================================================
# push.sh
# One-shot: commit + tag v1.0.0 + push + create GitHub release
# Run from repo root:
#   cd IBD-complement-analysis
#   bash push.sh
# Prereqs: gh CLI authenticated (`gh auth login`)
# =============================================================================
set -euo pipefail

TAG="v1.0.0"
REPO="zhangweiminshanghai/IBD-complement-analysis"
MSG="Supplementary Figures S1–S13: Gut-style overhaul (v1.0.0)"

echo "==> checking gh auth"
gh auth status || { echo "Run: gh auth login"; exit 1; }

echo "==> staging files"
git add CHANGELOG_Supplementary_Figures.md \
        code/R/*.R code/python/*.py \
        figures/*.png 2>/dev/null || true
git status --short

echo "==> committing"
git commit -m "$MSG" || echo "(nothing to commit)"

echo "==> tagging $TAG"
git tag -f "$TAG"
git push origin master --tags

echo "==> creating GitHub release"
gh release create "$TAG" \
    --repo "$REPO" \
    --title "Supplementary Figures v1.0.0 (Gut style)" \
    --notes "$(cat <<'EOF'
## Summary
Full overhaul of Supplementary Figures S1–S13 to meet *Gut* (BMJ) standards:
single-paragraph captions, British spelling, panel labels unified, figure dimensions ≤9×7 in, aspect-error <3 %, zero overlap (S10 v6 dual-axis, S11/S12/S13 rebuilt).

## Key fixes
| Fig | Change |
|---|---|
| S1 | New 6-tier vertical flowchart, core-narrative box "C2 initiated, CFI compensatory" |
| S2 | Caption merged to single Gut-style paragraph |
| S3 | Added panel D (PPI network of yellow-module hubs) |
| S6 | ROC curves with Hanley–McNeil 95 % CI bands |
| S10 | Dual-axis ridge plot — zero label/curve overlap |
| S11 | Binding-energy bars x-max 12, de-overlapped labels |
| S12 | Proper Radial MR placeholder, titles pixel-aligned |
| S13 | C-index panel relabelled, single-paragraph caption |

## Repo layout
\`\`\`
CHANGELOG_Supplementary_Figures.md   ← full modification report
code/R/                               ← R implementations
code/python/                           ← Python equivalents
figures/                               ← rendered PNGs (220 DPI)
push.sh                                ← this script
\`\`\`
EOF
)" \
    figures/S1_flowchart.png \
    figures/S2_batch_qc.png \
    figures/S3_wgcna.png \
    figures/S6_roc.png \
    figures/S10_gsea_v6.png \
    figures/S11_docking.png \
    figures/S12_mr.png \
    figures/S13_qc.png \
    CHANGELOG_Supplementary_Figures.md

echo "==> done. Release URL:"
gh release view "$TAG" --repo "$REPO" --json url -q .url
