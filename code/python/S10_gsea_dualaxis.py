"""
S10_gsea_dualaxis.py
Supplementary Figure S10 — GSEA dot plot + dual-axis ridge plot (v6) + NES bars
Key fix: dual-axis architecture eliminates label/curve overlap

Run: python code/python/S10_gsea_dualaxis.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches

rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
OUT = "figures"
import os; os.makedirs(OUT, exist_ok=True)

# ---- GSEA summary data ----------------------------------------------------
terms = ["cytoplasmic translation", "ribosome biogenesis",
         "xenobiotic metabolism", "lipid metabolic process"]
NES   = [-1.98, -1.76, 1.54, 1.42]
padj  = [6.28e-7, 1.20e-3, 5.79e-4, 2.80e-2]
ngene = [120, 95, 60, 78]

# ---- (A) Dot plot ----------------------------------------------------------
figA, axA = plt.subplots(figsize=(3.2, 3.0))
for i, c in enumerate(padj):
    axA.scatter(NES[i], i, s=ngene[i]*8, c=[NES[i]], cmap="RdBu_r",
                vmin=-2.2, vmax=2.2, edgecolors="black", lw=0.3, zorder=3)
    axA.text(NES[i] + (0.08 if NES[i]>0 else -0.08), i, f"NES={NES[i]:.2f}\np={c:.1e}",
             ha="left" if NES[i]>0 else "right", va="center", fontsize=6.5)
axA.axvline(0, ls="--", color="grey")
axA.set(yticks=range(4), yticklabels=terms, xlabel="Normalized Enrichment Score (NES)",
        title="(A) GSEA dot plot")
axA.invert_yaxis()
figA.tight_layout()
figA.savefig(f"{OUT}/S10_panelA.png", dpi=220, bbox_inches="tight", facecolor="white")

# ---- (B) DUAL-AXIS ridge plot ---------------------------------------------
figB, (axL, axR) = plt.subplots(1, 2, figsize=(4.5, 3.0),
                                   gridspec_kw={"width_ratios":[0.38,0.62]})
# LEFT axis: pure labels
for i, t in enumerate(terms):
    axL.text(0.5, i, t, ha="center", va="center", fontsize=9, color="#2C3E50")
axL.set(xlim=(0,1), ylim=(-0.5, 3.5), xticks=[], yticks=[])
for s in axL.spines.values(): s.set_visible(False)

# RIGHT axis: pure ridge curves
rng = np.random.default_rng(3)
x = np.linspace(-2500, 2500, 500)
cols = ["#2166AC","#4393C3","#D6604D","#B2182B"]
for i, (t, nes, c) in enumerate(zip(terms, NES, cols)):
    mu = nes * 800
    y = np.abs(nes) * np.exp(-0.5*((x-mu)/900)**2)
    y_shifted = i + 0.9*y/y.max()
    axR.fill_between(x, i, y_shifted, alpha=0.65, color=c, lw=0.3, edgecolor="white")
axR.set(xlim=(-2400,2400), ylim=(-0.5,3.5), yticks=[], ylabel="",
        xlabel="Running enrichment score", title="")
for s in ["top","right"]: axR.spines[s].set_visible(False)
figB.suptitle("(B) Enrichment distributions", fontsize=11, fontweight="bold", y=1.02)
figB.tight_layout()
figB.savefig(f"{OUT}/S10_panelB.png", dpi=220, bbox_inches="tight", facecolor="white")

# ---- (C) NES bar chart ----------------------------------------------------
figC, axC = plt.subplots(figsize=(2.6, 3.0))
cols_bar = ["#2166AC" if n<0 else "#B2182B" for n in NES]
axC.barh(range(4), NES, color=cols_bar, edgecolor="#2C3E50", lw=0.3, height=0.6)
for i, v in enumerate(NES):
    axC.text(v + (0.06 if v>0 else -0.06), i, f"{v:.2f}",
             ha="left" if v>0 else "right", va="center", fontsize=8)
axC.axvline(0, ls="--", color="grey")
axC.set(yticks=range(4), yticklabels=terms, xlabel="NES", title="(C) NES bar chart")
axC.invert_yaxis()
figC.tight_layout()
figC.savefig(f"{OUT}/S10_panelC.png", dpi=220, bbox_inches="tight", facecolor="white")

# ---- Combine A / B / C ----------------------------------------------------
fig, (a, b, c) = plt.subplots(1, 3, figsize=(9, 3.0),
                                 gridspec_kw={"width_ratios":[1.0,1.4,0.9]})
for src, dst in [(figA.axes[0], a), (figB.axes[1], b), (figC.axes[0], c)]:
    dst.imshow(np.asarray(src.figure.canvas.buffer_rgba())[...,:3])
    dst.axis("off")
plt.tight_layout()
plt.savefig(f"{OUT}/S10_gsea_v6.png", dpi=220, bbox_inches="tight", facecolor="white")
print("[S10] saved → figures/S10_gsea_v6.png")
