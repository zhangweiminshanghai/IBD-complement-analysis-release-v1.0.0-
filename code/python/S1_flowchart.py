"""
S1_flowchart.py
Supplementary Figure S1 — Vertical 6-tier flowchart (L1→L6)
Core narrative box: "C2 initiated, CFI compensatory"

Run: python code/python/S1_flowchart.py
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams

rcParams.update({"font.family":"DejaVu Sans","font.size":9})
OUT="figures"; import os; os.makedirs(OUT, exist_ok=True)

fig, ax = plt.subplots(figsize=(9,7.0))
ax.set_xlim(-6,6); ax.set_ylim(-1.2,9.2)
ax.axis("off")

# ---- tier definitions ----------------------------------------------------
tiers = [
  ("L1", 8.0, "#D6EAF8", "#2980B9",
   "Data acquisition & preprocessing\n• GSE16879 (discovery, n=85)\n• GSE75214 (validation, n=97)\n• Normalisation & batch correction (limma)"),
  ("L2a", 6.6, "#FADBD8", "#C0392B",
   "Differential expression\n• limma-voom\n• |log2FC|>1, FDR<0.05\n• Up/down genes → pathway analysis"),
  ("L2b", 5.2, "#FADBD8", "#C0392B",
   "WGCNA co-expression\n• Soft-threshold β=14\n• Yellow module\n• CFI, S100A8, S100A9, PAQR5, KCNE3"),
  ("L3", 3.8, "#E8DAEF", "#8E44AD",
   "Mendelian randomisation\n• C2 (protective) → CD\n• IVW OR=0.43 (95% CI 0.30–0.61)\n• p=2×10⁻⁶"),
  ("L4", 2.4, "#D5F5E3", "#27AE60",
   "Single-cell RNA-seq\n• GSE134809 (Seurat)\n• 160,981 cells, 12 cell types\n• Cell-type-specific complement"),
  ("L5a", 1.0, "#FCF3CF", "#D4AC0D",
   "3-gene signature\n• LASSO → 10-fold CV\n• CFI · PAQR5 · KCNE3\n• AUC=0.90 (disc.); AUC=0.98 (val.)"),
  ("L5b", -0.3, "#FCF3CF", "#D4AC0D",
   "Structural pharmacology\n• Nafamostat mesylate docking\n• C1r/C1s vs CFI selectivity\n• Therapeutic implication"),
]

def draw_tier(y, fill, border, txt):
    ax.add_patch(mpatches.FancyBboxPatch((-5,y-0.55),10,1.1,
                  boxstyle="round,pad=0.15",fc=fill,ec=border,lw=1.8))
    ax.text(0,y,txt,ha="center",va="center",fontsize=8.5,color="#2C3E50",
            linespacing=1.4)

for name, y, fc, ec, txt in tiers:
    draw_tier(y, fc, ec, txt)
    ax.text(-5.6,y,name,ha="right",va="center",fontsize=9,fontweight="bold",color="#2C3E50")

# arrows
for y1,y2 in [(7.25,6.95),(5.85,5.55),(4.45,4.15),(3.05,2.75),(1.55,1.25),(0.25,-0.05)]:
    ax.annotate("",xy=(0,y2+0.55),xytext=(0,y1-0.55),
                arrowprops=dict(arrowstyle="-|>",lw=1.3,color="#7F8C8D"))

# ---- core narrative box ---------------------------------------------------
ax.add_patch(mpatches.FancyBboxPatch((2.2,3.65),3.4,0.95,
              boxstyle="round,pad=0.2",fc="#C0392B",ec="#922B21",lw=2.2))
ax.text(3.9,4.12,"Core narrative:",ha="center",va="center",
        fontsize=10,fontweight="bold",color="white")
ax.text(3.9,3.80,"C2 initiated, CFI compensatory",ha="center",va="center",
        fontsize=9,fontstyle="italic",color="#FADBD8")

# ---- legend colour bar ----------------------------------------------------
legend_items = [
  ("Data / QC","#D6EAF8","#2980B9"),
  ("DE / WGCNA","#FADBD8","#C0392B"),
  ("MR","#E8DAEF","#8E44AD"),
  ("scRNA-seq","#D5F5E3","#27AE60"),
  ("Signature / Docking","#FCF3CF","#D4AC0D"),
]
for i,(nm,fc,ec) in enumerate(legend_items):
  ax.add_patch(mpatches.FancyBboxPatch((-5.8+i*1.3,-1.05),1.1,0.35,
                boxstyle="round,pad=0.08",fc=fc,ec=ec,lw=1.0))
  ax.text(-5.8+i*1.3+0.55,-0.87,nm,ha="center",va="center",fontsize=6.5,color="#2C3E50")

ax.set_title("(A) Multi-omics analytical pipeline", fontsize=12, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT}/S1_flowchart.png", dpi=220, bbox_inches="tight", facecolor="white")
print("[S1] saved → figures/S1_flowchart.png")
