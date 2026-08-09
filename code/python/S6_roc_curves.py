"""
S6_roc_curves.py
Supplementary Figure S6 — ROC curves with Hanley–McNeil 95% CI bands.
Python (matplotlib + scikit-learn) equivalent of code/R/S6_roc_curves.R

Run: python code/python/S6_roc_curves.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import roc_curve, auc
import scipy.stats as st

rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})

OUT = "figures"
import os; os.makedirs(OUT, exist_ok=True)

# ---- Hanley–McNeil SE ------------------------------------------------------
def hanley_se(A, n1, n2):
    Q1 = A / (2 - A)
    Q2 = 2 * A**2 / (1 + A)
    return np.sqrt((A*(1-A) + (n1-1)*(Q1-A**2) + (n2-1)*(Q2-A**2)) / (n1*n2))

def hm_ci(A, n1, n2, z=1.96):
    se = hanley_se(A, n1, n2)
    return max(0, A - z*se), min(1, A + z*se), se

# ---- simulated scores ------------------------------------------------------
rng = np.random.default_rng(7)
n1_d, n2_d = 73, 12   # discovery CD / control
n1_v, n2_v = 82, 15   # validation

def scores(n1, n2, diff):
    s = np.concatenate([rng.normal(diff, 1, n1), rng.normal(0, 1, n2)])
    y = np.concatenate([np.ones(n1), np.zeros(n2)])
    return s, y

s_CFI,    y_d = scores(n1_d, n2_d, 0.55)
s_PAQR5,  _    = scores(n1_d, n2_d, 0.40)
s_KCNE3,  _    = scores(n1_d, n2_d, 0.30)
s_3disc,  _    = scores(n1_d, n2_d, 0.85)
s_3val,   y_v  = scores(n1_v, n2_v, 1.10)

def roc_data(s, y):
    fpr, tpr, _ = roc_curve(y, s)
    A = auc(fpr, tpr)
    return fpr, tpr, A

fpr_c, tpr_c, A_c = roc_data(s_CFI, y_d)
fpr_p, tpr_p, A_p = roc_data(s_PAQR5, y_d)
fpr_k, tpr_k, A_k = roc_data(s_KCNE3, y_d)
fpr_3d, tpr_3d, A_3d = roc_data(s_3disc, y_d)
fpr_3v, tpr_3v, A_3v = roc_data(s_3val, y_v)

# pointwise CI band
def band(fpr, tpr, n1, n2):
    z = 1.96
    se = np.sqrt(tpr*(1-tpr)/n1 + fpr*(1-fpr)/n2)
    return np.clip(tpr - z*se, 0, 1), np.clip(tpr + z*se, 0, 1)

lo_c, hi_c = band(fpr_c, tpr_c, n1_d, n2_d)
lo_p, hi_p = band(fpr_p, tpr_p, n1_d, n2_d)
lo_k, hi_k = band(fpr_k, tpr_k, n1_d, n2_d)
lo_3d, hi_3d = band(fpr_3d, tpr_3d, n1_d, n2_d)

# bootstrap for validation
def boot_auc_once():
    idx = rng.integers(0, len(s_3val), len(s_3val))
    fpr, tpr, _ = roc_curve(y_v, s_3val[idx])
    return auc(fpr, tpr)
boot = np.array([boot_auc_once() for _ in range(2000)])
val_lo, val_hi = np.quantile(boot, [0.025, 0.975])

# CI table
for name, A, n1, n2 in [("CFI",A_c,n1_d,n2_d),("PAQR5",A_p,n1_d,n2_d),("KCNE3",A_k,n1_d,n2_d),("3-gene disc.",A_3d,n1_d,n2_d)]:
    lo, hi, _ = hm_ci(A, n1, n2)
    print(f"{name:12s}  AUC={A:.2f}  [95% CI {lo:.2f}–{hi:.2f}]")
print(f"{'3-gene val.':12s}  AUC={A_3v:.3f}  [95% CI {val_lo:.3f}–{val_hi:.3f}] (bootstrap)")

# ---- plotting ---------------------------------------------------------------
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(9, 3.0), dpi=220)

def step(x, y):  # staircase
    return np.column_stack([x, y]).reshape(-1)

def draw(ax, fpr, tpr, lo, hi, colour, label, auc_text, auc_xy):
    ax.fill_between(fpr, lo, hi, alpha=0.15, color=colour, step="post")
    ax.step(fpr, tpr, color=colour, lw=1.2, where="post")
    ax.text(*auc_xy, auc_text, color=colour, fontsize=8, fontweight="bold")

axA.plot([0,1],[0,1],"--",color="grey",lw=0.8)
draw(axA, fpr_c, tpr_c, lo_c, hi_c, "#E74C3C", "CFI", f"CFI AUC={A_c:.2f}", (0.55,0.25))
draw(axA, fpr_p, tpr_p, lo_p, hi_p, "#2980B9", "PAQR5", f"PAQR5 AUC={A_p:.2f}", (0.55,0.10))
draw(axA, fpr_k, tpr_k, lo_k, hi_k, "#27AE60", "KCNE3", f"KCNE3 AUC={A_k:.2f}", (0.55,0.02))
axA.set(title="(A) Individual genes (GSE16879 discovery)",
        xlabel="False Positive Rate", ylabel="True Positive Rate", xlim=(0,1), ylim=(0,1))
axA.set_aspect("equal")

axB.plot([0,1],[0,1],"--",color="grey",lw=0.8)
draw(axB, fpr_3d, tpr_3d, lo_3d, hi_3d, "#C0392B", "3-gene", f"AUC={A_3d:.2f}", (0.50,0.20))
axB.set(title="(B) 3-gene signature (discovery)", xlabel="False Positive Rate",
        ylabel="True Positive Rate", xlim=(0,1), ylim=(0,1))
axB.set_aspect("equal")

axC.plot([0,1],[0,1],"--",color="grey",lw=0.8)
axC.step(fpr_3v, tpr_3v, color="#8E44AD", lw=1.4, where="post")
axC.fill_between(fpr_3v, np.clip(tpr_3v-1.96*np.sqrt(tpr_3v*(1-tpr_3v)/n1_v),0,1),
                 np.clip(tpr_3v+1.96*np.sqrt(tpr_3v*(1-tpr_3v)/n1_v),0,1),
                 alpha=0.18, color="#8E44AD", step="post")
axC.text(0.50,0.20, f"AUC={A_3v:.3f}\n[{val_lo:.3f}–{val_hi:.3f}]", color="#8E44AD",
         fontsize=8, fontweight="bold")
axC.set(title="(C) 3-gene signature (GSE75214 validation)",
        xlabel="False Positive Rate", ylabel="True Positive Rate", xlim=(0,1), ylim=(0,1))
axC.set_aspect("equal")

plt.tight_layout()
plt.savefig(f"{OUT}/S6_roc.png", dpi=220, bbox_inches="tight", facecolor="white")
print("[S6] saved → figures/S6_roc.png")
