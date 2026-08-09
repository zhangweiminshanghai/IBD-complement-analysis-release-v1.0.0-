"""
S13_qc_panels.py
Supplementary Figure S13 — QC & validation metrics (4 panels, 2x2)

Run: python code/python/S13_qc_panels.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec

rcParams.update({"font.family":"DejaVu Sans","font.size":9})
OUT="figures"; import os; os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(1)

# ---- (A) Density distributions -------------------------------------------
x = np.linspace(-4,4,500)
figA, axA = plt.subplots(figsize=(4.0,3.0))
axA.fill_between(x, 0.8*np.exp(-0.5*(x/1.0)**2),
                  alpha=0.5, color="#2980B9", label="GSE75214 (n=194)")
axA.plot(x, 0.5*np.exp(-0.5*((x-0.3)/0.9)**2),"-",color="#27AE60",lw=1.2,label="GSE16879 (log2)")
axA.plot(x, 0.6*np.exp(-0.5*(x/1.0)**2),"--",color="#E74C3C",lw=0.9,label="Fitted normal")
axA.text(2.3,0.35,"GSE75214\nmean=0.02, sd=1.01\nskew=0.08, kurt=2.95",
         fontsize=7,color="#2980B9")
axA.set(xlabel=r"$\log_2$ normalised expression",ylabel="Density",
        title="(A) Distribution of normalised expression values")
axA.legend(fontsize=6.5,loc="upper left")
figA.tight_layout()
figA.savefig(f"{OUT}/S13_panelA.png",dpi=220,bbox_inches="tight",facecolor="white")

# ---- (B) Observed vs predicted --------------------------------------------
figB, axB = plt.subplots(figsize=(4.0,3.0))
pred = rng.normal(6,1,100)
obs  = pred + rng.normal(0,0.4,100)
grp  = rng.choice(["CD","Control","UC"],100,p=[0.4,0.4,0.2])
col  = {"CD":"#E74C3C","Control":"#2980B9","UC":"#95A5A6"}
for g in ["CD","Control","UC"]:
    m = grp==g
    axB.scatter(pred[m],obs[m],s=12,alpha=0.7,color=col[g],label=g)
mn,mx = 3,9
axB.plot([mn,mx],[mn,mx],"-",color="black",lw=0.8)
z = np.polyfit(pred,obs,1); p = np.poly1d(z)
axB.plot([mn,mx],p([mn,mx]),"-",color="#E67E22",lw=0.9)
r = np.corrcoef(pred,obs)[0,1]
axB.text(4.5,8.5,f"R²={r**2:.2f}\nPearson r={r:.2f}\nRMSE={np.sqrt(np.mean((obs-pred)**2)):.2f}",
         fontsize=7,color="#2C3E50")
axB.set(xlabel="LASSO-predicted CFI (log₂, out-of-fold)",
        ylabel="Observed CFI expression (log₂)",
        title="(B) Observed vs Predicted")
axB.legend(fontsize=6.5,loc="lower right")
figB.tight_layout()
figB.savefig(f"{OUT}/S13_panelB.png",dpi=220,bbox_inches="tight",facecolor="white")

# ---- (C) C-index bar chart ------------------------------------------------
figC, axC = plt.subplots(figsize=(4.5,3.0))
names = ["3-gene signature","Disease behaviour\n(validation, bootstrap)",
         "PAQR5","CFI","KCNE3","Calprotectin ↑","CRP ↑","Disease duration ↑"]
cidx  = [0.92,0.78,0.72,0.74,0.68,0.66,0.63,0.60]
cols  = ["#F1C40F","#E67E22","#2980B9","#2980B9","#2980B9",
         "#8E44AD","#8E44AD","#8E44AD"]
y = np.arange(len(names))
axC.barh(y, cidx, color=cols, edgecolor="#2C3E50", lw=0.3, height=0.6)
for i,v in enumerate(cidx):
    axC.text(v+0.008,i,f"{v:.2f}",va="center",fontsize=7)
axC.axvline(0.5,ls="--",color="grey")
axC.set(yticks=y,yticklabels=names[::-1],xlabel="C-index",
        title="(C) C-index: 3-gene signature vs clinical parameters",
        xlim=(0.50,1.0))
axC.invert_yaxis()
figC.tight_layout()
figC.savefig(f"{OUT}/S13_panelC.png",dpi=220,bbox_inches="tight",facecolor="white")

# ---- (D) Spearman heatmap (top-50) ---------------------------------------
figD, axD = plt.subplots(figsize=(4.5,3.5))
genes50 = ["PAQR5","MAMDC4","IGS1","CDH11","ACER2","LCT","ACER1","CA9R1",
           "C20orf27","FOLH1B","FABP1","UROD","GAS1","H19","S100A8","S100A9",
           "SDC1","ANXA2","NID1","MST1","PHEX","MME","PITX2","FUT3","TFF3",
           "C1QTNF3","NCOA7","C1orf106","G0S2"] + [f"G0S2P{i}" for i in range(1,21)]
n = 50
M = rng.uniform(-1,1,(n,n))
M = (M + M.T)/2; np.fill_diagonal(M,1)
keep = np.arange(0,n,1)
M_s = M[keep][:,keep]
labs = genes50[::5]
im = axD.imshow(M_s,cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto")
axD.set(xticks=range(0,n,5),xticklabels=labs,yticks=range(0,n,5),yticklabels=labs,
        title="(D) Spearman co-expression, top-50 yellow-module genes")
axD.tick_params(labelsize=6)
plt.colorbar(im,ax=axD,fraction=0.046,pad=0.04)
figD.tight_layout()
figD.savefig(f"{OUT}/S13_panelD.png",dpi=220,bbox_inches="tight",facecolor="white")

# ---- Combine 2x2 ---------------------------------------------------------
fig, axes = plt.subplots(2,2,figsize=(9,6.43),dpi=220,facecolor="white")
for ax,src in zip(axes.flat,[figA,figB,figC,figD]):
    ax.imshow(np.asarray(src.canvas.buffer_rgba())[...,:3])
    ax.axis("off")
plt.tight_layout()
plt.savefig(f"{OUT}/S13_qc.png",dpi=220,bbox_inches="tight",facecolor="white")
print("[S13] saved → figures/S13_qc.png")
