# IBD-complement-analysis

**C2 protects against Crohn's disease while CFI marks ileal-specific compensation: a multi-omics and Mendelian randomisation study**

This repository contains the analysis code, processed data, and publication-ready
figures for the manuscript above. It is released to satisfy the reproducibility
requirements of the *Gut* review: all scripts, the three previously-missing
analysis scripts, environment specifications, and 300-DPI vector figures are
provided so that every result in the paper can be regenerated.

> **Double-blind review.** This is the anonymous review copy. The corresponding
> author's identity, affiliation, and any identifying metadata have been removed
> from code comments. **Do not share the raw GitHub URL** — it exposes the owner
> handle and commit-author email.
>
> ### ✅ Anonymous link for reviewers (use this)
> **https://anonymous.4open.science/r/IBD-complement-analysis-release-v1.0.0--7783/**
>
> This mirror hides the GitHub owner, commit-author name/email, and original
> repository path. Submit **only this URL** to the journal.
>
> After acceptance, the versioned release (`v1.0.0`) and GitHub handle can be
> revealed publicly.

---

## 1. Repository structure

```
IBD-complement-analysis/
├── README.md                     # this file
├── LICENSE                       # MIT
├── .gitattributes                # Git LFS rules for figures / matrices
├── environment/
│   ├── renv.lock                 # R package versions (R 4.5.3)
│   └── conda.yml                 # conda env (R + Python + AutoDock Vina)
├── data/
│   ├── README.md                 # GEO accessions + processing notes
│   ├── extract_phenotype.py      # regenerate phenotype CSVs from GEO matrix
│   └── processed/                # processed matrices + phenotypes (LFS)
├── scripts/
│   ├── 01_DEG_limma.R           # limma DEG (GSE16879 / GSE75214)
│   ├── 02_WGCNA.R               # WGCNA co-expression modules
│   ├── 03_LASSO_RF.R            # LASSO + Random Forest 3-gene signature
│   ├── 04_CIBERSORT.R           # CIBERSORT immune deconvolution
│   ├── 05_scRNA_Seurat.R        # Seurat v5 single-cell CFI+ map
│   ├── 06_MR_TwoSampleMR.R       # TwoSampleMR (C2 protective, CFI null)
│   ├── 07_docking_AutoDockVina.py# AutoDock Vina wrapper (nafamostat)
│   ├── 08_bootstrap_AUC.R       # bootstrap signature AUC (2000 resamples)
│   ├── CFI_stratified_DEG.R     # ★ median-split CFI DEG (CFI-pathway)
│   ├── three_tier_CFI_C2.R       # ★ CFI–C2 three-tier framework
│   ├── docking_postprocess.py    # ★ parse Vina output -> DeltaG + contacts
│   └── manuscript.Rmd           # manuscript RMarkdown (knit-safe)
├── results/
│   └── figures/                  # every submitted figure: .png(300dpi)+.pdf (+.svg)
└── docs/
    ├── ST1_batch_effects.pdf     # batch-effect PCA (before/after correction)
    ├── ST2_WGCNA_modules.pdf     # WGCNA dendrogram + module-trait heatmap
    └── Figure_S1_workflow.pdf    # integrative multi-omics workflow
```
★ = scripts that were missing from the original submission and are now supplied.

---

## 2. Software versions

| Component | Version |
|-----------|---------|
| R | 4.5.3 |
| Bioconductor | 3.20 |
| Seurat | 5.5.0 |
| limma | 3.66.0 |
| WGCNA | 1.74 |
| glmnet | 5.0 |
| randomForest | 4.7.1.2 |
| TwoSampleMR | 0.7.4 |
| clusterProfiler | 4.18.4 |
| Python | 3.10+ |
| AutoDock Vina | 1.2.3 (optional; for Figure 6 regeneration) |
| Git LFS | required to clone figure/matrix binaries |

**Random seeds:** all stochastic steps use `set.seed(2024)` (R) /
`random_state=2024` (Python) for exact reproducibility.

---

## 3. Environment setup

**R users** (recommended):
```r
install.packages("renv")
renv::restore()                 # uses environment/renv.lock
```

**Conda users:**
```bash
conda env create -f environment/conda.yml
conda activate ibd-complement-analysis
```

---

## 4. Data

Processed matrices only are committed (`data/processed/`). Raw data are retrieved
from GEO (no raw FASTQ/CEL in this repo):

| Cohort | GSE | Platform | Role |
|--------|-----|----------|------|
| Discovery | GSE16879 | GPL570 | DEG / signature training (73 CD, 12 control) |
| Validation | GSE75214 | GPL11532 | DEG / signature validation (75 CD, 22 control) |
| scRNA-seq | GSE134809 | 10x Genomics | CFI+ single-cell map (160,981 cells) |

Phenotype CSVs can be regenerated from the GEO series matrix:
```bash
python data/extract_phenotype.py
```

---

## 5. Script ↔ figure / table mapping

| Script | Produces | Manuscript element |
|--------|----------|--------------------|
| `01_DEG_limma.R` | volcano plots, DEG tables | Figure 1A–B, Table 2 |
| `02_WGCNA.R` | module dendrogram, trait heatmap | `docs/ST2_WGCNA_modules.pdf` |
| `03_LASSO_RF.R` | LASSO/RF, ROC | Figure 2 |
| `04_CIBERSORT.R` | immune fractions | Figure 3 |
| `05_scRNA_Seurat.R` | UMAP, CFI+ proportions | Figure 4 |
| `06_MR_TwoSampleMR.R` | C2 OR=0.43, CFI p=0.73 | Results (MR) |
| `07_docking_AutoDockVina.py` + `docking_postprocess.py` | ΔG, contacts | Figure 6 |
| `08_bootstrap_AUC.R` | AUC 0.711 (95% CI 0.565–0.844) | Table 2 |
| `CFI_stratified_DEG.R` | CFI-high vs low DEG | Figure 5 (CFI-pathway) |
| `three_tier_CFI_C2.R` | 3-tier framework table | Discussion |

---

## 6. Reproduce the key figures

```r
# Figure 2 - machine-learning signature (LASSO + RF + ROC)
Rscript scripts/03_LASSO_RF.R

# Figure 4 - single-cell CFI+ map (Seurat v5; needs GSE134809)
Rscript scripts/05_scRNA_Seurat.R

# Figure 6 - nafamostat docking (needs AutoDock Vina + 2XRC_clean.pdb)
python scripts/07_docking_AutoDockVina.py --receptor 2XRC_clean.pdb \
       --ligand nafamostat.sdf --outdir results/figures
python scripts/docking_postprocess.py results/figures/docking_out.pdbqt
```

All other figures are regenerated by running `01` → `08` in order; outputs land
in `results/figures/` as 300-DPI PNG + vector PDF. The committed figures in
`results/figures/` are the exact versions used in the submission.

---

## 7. Key results recap

- **C2** is genetically protective against CD (MR OR = 0.43 per SD, p = 2.2×10⁻⁶).
- **CFI** shows no independent causal effect (p = 0.73) → a downstream
  compensatory marker, enriched in ileal CD single cells (0.89%) and absent in
  colonic CD.
- A **3-gene signature** (CFI, PAQR5, KCNE3) discriminates CD with
  cross-validated AUC = 0.90; bootstrap AUC = 0.711 (95% CI 0.565–0.844).
- **Structural docking** shows nafamostat binds C1r/C1s (ΔG ≈ −8.8 kcal/mol)
  with ~2 kcal/mol selectivity over CFI (ΔG ≈ −6.8), repositioning nafamostat
  as a structural probe rather than a CD therapy.

---

## 8. License & citation

Code is released under the MIT License (see `LICENSE`). If you use this work,
please cite the manuscript. Data are redistributed under the original GEO
accession terms.
