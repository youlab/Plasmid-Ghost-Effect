import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from calculate_growth_rate import GR

colors = ["#AAAAAA", "#B53030", "#FF9000"]

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
    "DA28102": "./processed_data/DA28102_ancestor_LB+Cm_OD.xlsx",
    "pCU1": "./processed_data/pCU1_ancestor_LB+Cm_OD.xlsx",
    "R6K": "./processed_data/R6K_ancestor_LB+Cm_OD.xlsx",
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

# growth-rate calculation from saved OD xlsx (one value per replicate)
def growth_rate_from_replicates(time: np.ndarray, od_matrix: np.ndarray, labels: list[str], plasmid: str) -> pd.DataFrame:
    records = []
    for i, rep_label in enumerate(labels):
        match = re.fullmatch(r"c(\d+)_(\d+)", rep_label)
        if match is None:
            raise ValueError(f"Invalid replicate label '{rep_label}'. Expected format: cx_y")
        colony = int(match.group(1))
        biological_replicate = int(match.group(2))

        y_od = np.maximum(od_matrix[:, i], 0.05)
        logy_od = np.log(y_od)
        _, slope = GR(time, logy_od, window_size=12)
        records.append(
            {
                "plasmid": plasmid,
                "replicate": rep_label,
                "colony": colony,
                "biological_replicate": biological_replicate,
                "growth_rate": float(np.max(slope)),
            }
        )
    return pd.DataFrame(records)


df_gr = pd.concat(
    [
        growth_rate_from_replicates(time_da, da_od, da_reps, "DA28102"),
        growth_rate_from_replicates(time_pcu1, pcu1_od, pcu1_reps, "pCU1"),
        growth_rate_from_replicates(time_r6k, r6k_od, r6k_reps, "R6K"),
    ],
    ignore_index=True,
)
df_gr.to_excel("./processed_data/conjugative_growth_rates_ancestor.xlsx", index=False)

# mean +/- sem summary bar plot
order = ["DA28102", "pCU1", "R6K"]

summary_rows = []
for plasmid in order:
    vals = df_gr.loc[df_gr["plasmid"] == plasmid, "growth_rate"].to_numpy(dtype=float)
    mean_val = float(np.mean(vals))
    sem_val = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
    summary_rows.append(
        {
            "plasmid": plasmid,
            "mean_growth_rate": mean_val,
            "sem_growth_rate": sem_val,
        }
    )
summary = pd.DataFrame(summary_rows)

# plasmid burden = max growth rate(plasmid) / max growth rate(DA28102)
# Here we use mean max-growth-rate across replicates and propagate SEM for a ratio:
# sigma_r = r * sqrt((sigma_x/x)^2 + (sigma_y/y)^2)
summary = summary.set_index("plasmid")
da_mean = np.array(summary.loc["DA28102", "mean_growth_rate"], dtype=float).item()
da_sem = np.array(summary.loc["DA28102", "sem_growth_rate"], dtype=float).item()

summary["burden"] = summary["mean_growth_rate"] / da_mean

rel_var = (summary["sem_growth_rate"] / summary["mean_growth_rate"]) ** 2 + (da_sem / da_mean) ** 2
summary["burden_sem"] = summary["burden"] * np.sqrt(rel_var)
#summary.loc[summary.index == "DA28102", "burden_sem"] = 1

summary.to_excel("./processed_data/conjugative_plasmid_burden_ancestor.xlsx", index=False)

fig_gr, ax_gr = plt.subplots(figsize=(3, 3))
x = np.arange(len(order))

error_kw = {"elinewidth": 1.2, "capsize": 4, "capthick": 1.2}
ax_gr.bar(x,summary["mean_growth_rate"].to_numpy(dtype=float),
    yerr=summary["sem_growth_rate"].fillna(0).to_numpy(dtype=float),
    width=0.6, error_kw=error_kw, lw=1.2, facecolor="None", edgecolor="k", zorder=2)

# scatters represent mean + sem for each colony.
rng = np.random.default_rng(42)
jitter_scale = 0.08
for j, plasmid in enumerate(order):
    df_p = df_gr.loc[df_gr["plasmid"] == plasmid]
    colony_ids = np.sort(df_p["colony"].unique())

    colony_means = []
    colony_sems = []
    for colony_id in colony_ids:
        colony_vals = np.array(df_p.loc[df_p["colony"] == colony_id, "growth_rate"], dtype=float)
        colony_means.append(float(np.mean(colony_vals)))
        colony_sems.append(float(np.std(colony_vals, ddof=1) / np.sqrt(len(colony_vals))) if len(colony_vals) > 1 else 0.0)

    y_vals = np.array(colony_means, dtype=float)
    y_err = np.array(colony_sems, dtype=float)
    x_j = rng.normal(loc=x[j], scale=jitter_scale, size=len(y_vals))
    ax_gr.scatter(x_j, y_vals, s=60, linewidths=1, facecolor=colors[j], edgecolors="k", zorder=3)
    ax_gr.errorbar(x_j, y_vals, yerr=y_err, marker="o", markersize=0, elinewidth=1.2,
        linewidth=0, capsize=2, capthick=1.2, color="k", zorder=4)

ax_gr.set_xticks(x)
ax_gr.set_xticklabels(order,rotation=90)
ax_gr.set_yticks([0, 0.5, 1.0])
ax_gr.set_ylabel("growth rate (1/hour)")
fig_gr.tight_layout()
fig_gr.savefig("./figures/conjugative_growth_rate_ancestor.png", dpi=300)

plt.show()