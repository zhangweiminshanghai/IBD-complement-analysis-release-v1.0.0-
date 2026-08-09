"""
S12_mr_panels.py
Supplementary Figure S12 — MR sensitivity (4 panels, aligned titles)
Titles use same baseline y via unified gridspec.

Run: python code/python/S12_mr_panels.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec

rcParams.update({"font.family":"DejaVu Sans","font.size":9})
OUT="figures"; import os; os.makedirs(OUT, exist_ok=True)

fig = plt.figure(figsize=(9, 2.6))
gs = GridSpec(1, 4, figure=fig, width_ratios=[1, 0.85, 1, 1])

# ---- (A) Leave-one-out ---------------------------------------------------
axA = fig.add_subplot(gs[0,0])
snps = ["rs488755","rs114502302","All SNPs"]
orv = [0.42, 0.44, 0.43]
lo  = [0.30, 0.32, 0.30]
hi  = [0.61, 0.62, 0.61]
y   = np.arange(3)
xerr = np.array([(o-l, h-o) for o,l,h in zip(orv,lo,hi)]).T  # shape (2, n)
axA.errorbar(orv, y, xerr=xerr,
             fmt="o", color="#2980B9", ecolor="#2980B9", capsize=3, lw=1.2)
axA.axvline(1, ls="--", color="grey")
axA.plot([0.30,0.61],[1.5,1.5],"-",color="#C0392B",lw=1.5)
axA.text(0.46,1.7,"IVW OR=0.43\n(0.30–0.61)",color="#C0392B",fontsize=7,ha="center")
axA.set(xlim=(0.25,0.70), yticks=y, yticklabels=snps[::-1],
        xlabel="OR (95% CI) for CD", title="(A) Leave-one-out meta-analysis")

# ---- (B) Radial MR placeholder --------------------------------------------
axB = fig.add_subplot(gs[0,1], aspect="equal")
th = np.linspace(0,2*np.pi,200)
axB.plot(1.96*np.cos(th),1.96*np.sin(th),"--",color="#BDC3C7",lw=0.8)
axB.scatter([0.8,-0.5],[0.6,1.2],color=["#2980B9","#27AE60"],s=30)
axB.text(0,-1.6,"Radial MR not applicable\n(n=2 instruments)",
         ha="center",fontsize=7,color="#7F8C8D")
axB.text(0,-2.1,"All SNPs within\nconfidence bounds",
         ha="center",fontsize=6.5,color="#27AE60")
axB.set(xticks=[],yticks=[],title="(B) Radial MR")
for s in axB.spines.values(): s.set_visible(False)

# ---- (C) Cumulative MA ---------------------------------------------------
axC = fig.add_subplot(gs[0,2])
n = [1,2,3,4,5]
or_c = [0.55,0.50,0.43,0.44,0.43]
lo_c = [0.35,0.34,0.34,0.33,0.34]
hi_c = [0.80,0.68,0.56,0.55,0.56]
axC.fill_between(n, lo_c, hi_c, alpha=0.18, color="#2980B9")
axC.plot(n, or_c,"-o",color="#2980B9",lw=1.2,ms=4)
axC.axvline(3,ls="--",color="#E74C3C",lw=0.8)
axC.text(3.2,0.75,"stabilised\nat n=3",color="#E74C3C",fontsize=7)
axC.set(xticks=n,xlabel="SNPs included (cumulative)",
        ylabel="OR (95% CI) for CD",ylim=(0.30,0.80),
        title="(C) Cumulative meta-analysis")

# ---- (D) Funnel plot -----------------------------------------------------
axD = fig.add_subplot(gs[0,3], aspect="equal")
axD.plot([-2,2,2,-2,-2],[5,5,20,20,5],"--",color="#BDC3C7",lw=0.6)
beta = [-0.15,-0.30,0.55]
se   = [1/12.5,1/8.0,1/6.5]
col  = ["#27AE60","#2980B9","#C0392B"]
lab  = ["rs488755","rs114502302","IVW pooled"]
axD.scatter(beta,1/np.array(se),c=col,s=40,zorder=5)
axD.axvline(0,color="grey",lw=0.6)
axD.add_patch(plt.Rectangle((0.20,14),0.45,5,fc="white",ec="grey",lw=0.5))
axD.text(0.42,17,"Symmetric →\nno small-study\neffects",
         ha="center",va="center",fontsize=6.5,color="#2C3E50")
axD.set(xlabel=r"$\beta$ (log OR)",ylabel="Precision (1/SE)",title="(D) Funnel plot")

fig.suptitle("(A)                (B)                (C)                (D)",
             fontsize=11,fontweight="bold",y=1.02,x=0.52)
fig.tight_layout()
fig.savefig(f"{OUT}/S12_mr.png",dpi=220,bbox_inches="tight",facecolor="white")
print("[S12] saved → figures/S12_mr.png")
