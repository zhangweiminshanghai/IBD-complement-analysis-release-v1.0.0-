# data/

Processed, analysis-ready matrices only. **No raw sequencing data is committed**
(raw FASTQ / CEL are retrieved from GEO; see accessions below).

## GEO accessions

| Cohort | GSE | Platform | Role | CD | Control | UC |
|--------|-----|----------|------|----|---------|----|
| Discovery | GSE16879 | GPL570 (HG-U133_Plus_2) | DEG / signature training | 73 | 12 | 48 |
| Validation | GSE75214 | GPL11532 (HuGene-1_0-st) | DEG / signature validation | 75 | 22 | 97 |
| scRNA-seq | GSE134809 | 10x Genomics | CFI+ single-cell map | 160,981 cells | – | – |

## Files in `processed/`

- `GSE16879_expr_normalized.txt` — log2-normalized expression matrix (probes→genes,
  hgu133plus2SYMBOL annotation), genes × samples. Produced by `scripts/01_DEG_limma.R`
  (pre-processing stage of `GEO_CD.R`).
- `GSE75214_expr_normalized.txt` — log2-normalized expression matrix (hugene10sttranscriptcluster
  annotation), genes × samples.
- `GSE16879_phenotype.csv` / `GSE75214_phenotype.csv` — `sample,group` (CD / Control / UC).
  Regenerated from the GEO series matrix with `extract_phenotype.py`
  (group inferred from `!Sample_source_name_ch1`: CDc/CDi/CD → CD, UC → UC, control → Control).
- `GSE16879_DEG_results.csv` / `GSE75214_DEG_results.csv` — full limma top-tables
  (logFC, AveExpr, t, P.Value, adj.P.Val, Change).
- `MR_C2_CFI.csv` — two-sample MR summary (C2 genetically protective OR=0.43, p=2.2×10⁻⁶;
  CFI null, p=0.73). Consumed by `scripts/three_tier_CFI_C2.R` (Tier 3).

## Reproduce the processed matrices from raw GEO

```r
# in R 4.4+ with GEOquery, limma, annotate
source("scripts/01_DEG_limma.R")   # performs download + normalize + annotate
```

Phenotype CSVs can be regenerated without R:

```bash
python data/extract_phenotype.py
```

## Notes

- Disease labels follow the manuscript: CD = Crohn's disease (ileal + colonic),
  UC = ulcerative colitis. DEG / ML signatures are trained on CD vs Control only
  (UC excluded). scRNA CFI+ analysis is restricted to ileal CD.
- For the single-cell tier of `three_tier_CFI_C2.R`, provide
  `scRNA_cfi_c2.csv` (columns: `cfi_positive` 0/1, `c2_expr` numeric) if you wish
  to compute Tier-1 statistics; otherwise Tier 1 is skipped with a note.
