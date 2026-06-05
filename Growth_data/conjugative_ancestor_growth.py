import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

colors = ["#AAAAAA",  "#74B3EB", "#FF9000"]

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
    "DA28102": "./processed_data/ancestor/DA28102_OD.xlsx",
    "pCU1": "./processed_data/ancestor/pCU1_OD.xlsx",
    "R6K": "./processed_data/ancestor/R6K_OD.xlsx",
}

for strain, xlsx_path in od_files.items():
    if not os.path.exists(xlsx_path):
        raise FileNotFoundError(f"Missing OD xlsx for {strain}: {xlsx_path}")

time_da, da_od, da_mean, da_std, da_reps = load_od_curve(od_files["DA28102"])
time_pcu1, pcu1_od, pcu1_mean, pcu1_std, pcu1_reps = load_od_curve(od_files["pCU1"])
time_r6k, r6k_od, r6k_mean, r6k_std, r6k_reps = load_od_curve(od_files["R6K"])

# growth curve comparison across strains, averaged across individual colonies and replicates
fig, ax = plt.subplots(figsize=(3, 3))
ax.plot(time_da, da_mean, linewidth=1.5, c=colors[0], label="DA28102")
ax.fill_between(time_da, da_mean - da_std, da_mean + da_std, color=colors[0], alpha=0.25)

ax.plot(time_pcu1, pcu1_mean, linewidth=1.5, c=colors[1], label="pCU1")
ax.fill_between(time_pcu1, pcu1_mean - pcu1_std, pcu1_mean + pcu1_std, color=colors[1], alpha=0.25)

ax.plot(time_r6k, r6k_mean, linewidth=1.5, c=colors[2], label="R6K")
ax.fill_between(time_r6k, r6k_mean - r6k_std, r6k_mean + r6k_std, color=colors[2], alpha=0.25)

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
fig.savefig("./figures/conjugative_growth_curve_ancestor.png", dpi=300)
fig.savefig("./figures/conjugative_growth_curve_ancestor.svg", dpi=300)
plt.show()