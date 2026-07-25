import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

cmap = plt.get_cmap('Set2')
colors = ["#AAAAAA"]+[cmap(i) for i in np.linspace(0, 1, 8)]

processed_dir = "./processed_data"

def load_od_curve(xlsx_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Return time, replicate matrix, mean OD, std OD, and replicate labels."""
    df = pd.read_excel(xlsx_path)
    time = df["time"].to_numpy(dtype=float)
    replicate_cols = [c for c in df.columns if c != "time"]
    od_values = df[replicate_cols].to_numpy(dtype=float)
    mean_od = np.mean(od_values, axis=1)
    std_od = np.std(od_values, axis=1)
    return time, od_values, mean_od, std_od, replicate_cols


# Files exported by growth_curve_extraction.py
od_files = {
    "MG1655": "./processed_data/ancestor/MG1655_OD.xlsx",
    "pSC101": "./processed_data/ancestor/pSC101_OD.xlsx",
    "colE1": "./processed_data/ancestor/colE1_OD.xlsx",
    "pUC": "./processed_data/ancestor/pUC_OD.xlsx",
}

for strain, xlsx_path in od_files.items():
    if not os.path.exists(xlsx_path):
        raise FileNotFoundError(f"Missing OD xlsx for {strain}: {xlsx_path}")

time_mg, mg_od, mg_mean, mg_std, mg_reps = load_od_curve(od_files["MG1655"])
time_psc101, psc101_od, psc101_mean, psc101_std, psc101_reps = load_od_curve(od_files["pSC101"])
time_cole1, cole1_od, cole1_mean, cole1_std, cole1_reps = load_od_curve(od_files["colE1"])
time_puc, puc_od, puc_mean, puc_std, puc_reps = load_od_curve(od_files["pUC"])

# growth curve comparison across strains, averaged across individual colonies and replicates
fig, ax = plt.subplots(figsize=(3, 3))
ax.plot(time_mg, mg_mean, linewidth=1.5, c=colors[0], label="MG1655")
ax.fill_between(time_mg, mg_mean - mg_std, mg_mean + mg_std, color=colors[0], alpha=0.25)

ax.plot(time_psc101, psc101_mean, linewidth=1.5, c=colors[1], label="pSC101")
ax.fill_between(time_psc101, psc101_mean - psc101_std, psc101_mean + psc101_std, color=colors[1], alpha=0.25)

ax.plot(time_cole1, cole1_mean, linewidth=1.5, c=colors[2], label="colE1")
ax.fill_between(time_cole1, cole1_mean - cole1_std, cole1_mean + cole1_std, color=colors[2], alpha=0.25)

ax.plot(time_puc, puc_mean, linewidth=1.5, c=colors[3], label="pUC")
ax.fill_between(time_puc, puc_mean - puc_std, puc_mean + puc_std, color=colors[3], alpha=0.25)

ax.set_yticks([0.0, 0.5, 1.0])
ax.set_xticks([0, 12, 24])
ax.set_xlabel("time (hours)")
ax.set_ylabel("OD600")
ax.legend(
    loc="lower right",
    fontsize=12,
    handlelength=0.8,
    handletextpad=0.4,
    labelspacing=0.35,
    frameon=False,
)

fig.tight_layout()
# fig.savefig("./figures/GFP_growth_curve_ancestor.png", dpi=300)
# fig.savefig("./figures/GFP_growth_curve_ancestor.svg", dpi=300)
plt.show()