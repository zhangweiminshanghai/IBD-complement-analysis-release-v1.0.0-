"""
S11_docking.py
Supplementary Figure S11 — Nafamostat binding energies (C-panel fixed)
*Gut* style: weak/strong groups, x-axis 0–12, no overlap

Run: python code/python/S11_docking.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
OUT = "figures"; import os; os.makedirs(OUT, exist_ok=True)

targets = ["TNFSF10","SGS","TYRO3","TRAF2","TIRAP","TNKS4"]
dG = [1.2, 1.5, 8.9, 9.6, 9.7, 9.9]
group = ["weak","weak","strong","strong","strong","strong"]
cols  = {"weak":"#BDC3C7","TYRO3":"#8E44AD","TRAF2":"#2980B9",
         "TIRAP":"#27AE60","TNKS4":"#C0392B"}

fig, ax = plt.subplots(figsize=(9, 3.2))
y_pos = list(range(6))
for i, (t, v, g) in enumerate(zip(targets, dG, group)):
    c = cols[g] if g=="weak" else cols[t]
    ax.barh(i, v, color=c, edgecolor="#2C3E50", lw=0.3, height=0.55)
    ax.text(v + 0.25, i, f"{v:.1f}", va="center", fontsize=9, color="#2C3E50")

ax.axvline(7, ls="--", color="#E74C3C", lw=0.8)
ax.text(7.3, 5.6, "strong binders →", color="#E74C3C", fontsize=9)
ax.text(9.7, 2.0, "C1r/C1s related", color="#2980B9", fontsize=8, style="italic")
ax.text(9.7, 1.0, "Weak/non-specific", color="#7F8C8D", fontsize=8, style="italic")
ax.set(xlim=(0,12), xlabel=r"Binding free energy  $\Delta$G (kcal/mol)",
       yticks=y_pos, yticklabels=targets[::-1],
       title="(C) Nafamostat binding interactions")
ax.invert_yaxis()
plt.tight_layout()
plt.savefig(f"{OUT}/S11_docking.png", dpi=220, bbox_inches="tight", facecolor="white")
print("[S11] saved → figures/S11_docking.png")
