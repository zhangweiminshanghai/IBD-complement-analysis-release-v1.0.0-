# IBD Figures Extraction — Summary

**Date:** 2026-08-03
**Repo:** `D:\IBD_Project\IBD-complement-analysis`
**Output dir:** `results/figures/`

## What was done

1. **Located final manuscript files** (v57):
   - Main: `Gut_manuscript_FINAL_v57.docx` → 13 embedded images
   - Supplementary: `Gut_Supplementary_FINAL_v57.docx` → 13 embedded images

2. **Extracted images from docx files** (read-only unzip):
   - Main docx: images named `Fig1A–D`, `Fig2A–D`, `Fig3A–C`, `Fig6A` based on XML position analysis
     - Note: main docx had NO formal "Figure N." caption paragraphs (all Normal style); figure numbers inferred from drawing positions relative to Results section text
     - image1.png = standalone workflow schematic in Methods (not a main figure)
   - Supp docx: images correctly named `SuppFig1–13` via `Figure S(N).` caption regex

3. **PDF generation**: every PNG → 300 DPI PDF via PIL

4. **SVG generation** (Rscript at `C:\Program Files\R\R-4.5.3\bin\Rscript.exe`):
   - ✅ Fig3.svg, Fig4.svg, Fig5.svg
   - ✅ SuppFig9.svg, SuppFig10.svg, SuppFig12.svg
   - R scripts for these 6 figures were also copied into `scripts/`

5. **Standalone file copies**:
   - `Figure_S1_Workflow_v26.{png,pdf}` → `SuppFig1.{png,pdf}`
   - `Figure_S13_CellTalk_v26.{png,pdf}` → `SuppFig13.{png,pdf}`
   - `Figure_Graphical_Abstract_v2.{png,pdf}` → `Graphical_Abstract.{png,pdf}`

6. **MANIFEST.csv** written with columns: figure, png_path, pdf_path, svg_or_source, source_script, notes

## Output summary

| Category | Count | Notes |
|---|---|---|
| PNG files | 26 | Fig1A–D, Fig2A–D, Fig3A–C, Fig6A (main docx); SuppFig1–13 (supp docx + standalone); Graphical_Abstract |
| PDF files | 27 | 300 DPI, one per PNG above |
| SVG files | 6 | Fig3, Fig4, Fig5, SuppFig9, SuppFig10, SuppFig12 |
| MANIFEST.csv | 1 | 27 rows |
| **Total files** | **59** | **7.15 MB** |

## Notable findings / gaps

- **Fig4 and Fig5 are NOT embedded** in the main docx. Their R scripts (`Figure4_main.R`, `Figure5_main.R`) exist in `D:\IBD_Project/` and have been copied to `scripts/`; standalone SVGs were successfully generated from those scripts.
- **Fig3 panels (A/B/C)**: PNGs extracted from docx; SVG generated from `Figure3_main.R` covers all 3 panels.
- **SuppFig1 and SuppFig13**: PDFs copied from standalone workflow/docking files (higher quality than docx-embedded versions).
- **image1.png from main docx** (Methods workflow schematic) was discarded as it is not a numbered figure panel.

## Key technical notes

- Main docx XML analysis: all 13 images are `<w:drawing>` with `r:embed`, all 13 images are in-media, but zero paragraphs have a non-Normal `w:pStyle` and zero contain the literal text "Figure N." (the captions are inline in the figure panel itself, not in the paragraph text). Figure numbers were inferred by mapping drawing positions to the Results section text layout.
