import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from calculate_growth_rate import GR

colors = ["#74B3EB", "#FF9000"]

processed_dir = "./processed_data"

def load_od_curve(xlsx_path: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return time, replicate matrix, and replicate labels."""
    df = pd.read_excel(xlsx_path)
    time = df["time"].to_numpy(dtype=float)
    replicate_cols = [c for c in df.columns if c != "time"]
    od_values = df[replicate_cols].to_numpy(dtype=float)
    return time, od_values, replicate_cols

# AUC calculation from saved OD xlsx (one value per replicate)
def quantify_fitness(time: np.ndarray, od_matrix: np.ndarray, labels: list[str], plasmid: str) -> pd.DataFrame:
    records = []
    for colony in range(1,4,1):
        auc = []
        gr = []
        for tech_rep in range(1,4,1):
            rep_label = f"c{colony}_{tech_rep}"
            od = od_matrix[:, labels.index(rep_label)]
            auc.append(float(np.sum((od[:-1] + od[1:]) * np.diff(time) * 0.5)))
            od = np.maximum(od, 0.01)  # floor to avoid log(0) and small value fluctuations
            gr.append(np.max(GR(time[7:], np.log(od[7:]), window_size=6)[1]))  # max growth rate from rolling window, starting after 30 minutes
        records.append(
            {
                "plasmid": plasmid,
                "time": 0,
                "replicate": 1,
                "colony": colony,
                "mean auc": np.mean(auc),
                "sem auc": np.std(auc, ddof=1) / np.sqrt(len(auc)),
                "mean gr": np.mean(gr),
                "sem gr": np.std(gr, ddof=1) / np.sqrt(len(gr))
            }
        )
    return pd.DataFrame(records)

# Files exported by growth_curve_extraction.py
od_files = {
    "DA28102": "./processed_data/ancestor/DA28102_OD.xlsx",
    "pCU1": "./processed_data/ancestor/pCU1_OD.xlsx",
    "R6K": "./processed_data/ancestor/R6K_OD.xlsx",
}

time_da, da_od, da_reps = load_od_curve(od_files["DA28102"])
time_pcu, pcu_od, pcu_reps = load_od_curve(od_files["pCU1"])
time_r6k, r6k_od, r6k_reps = load_od_curve(od_files["R6K"])

df = pd.concat(
    [
        quantify_fitness(time_da, da_od, da_reps, "DA28102"),
        quantify_fitness(time_pcu, pcu_od, pcu_reps, "pCU1"),
        quantify_fitness(time_r6k, r6k_od, r6k_reps, "R6K"),
    ],
    ignore_index=True,
)

# Calculate mean GR and AUC for DA28102 across its three colonies
da_df = df.loc[df["plasmid"] == "DA28102"]
da_colony_grs = da_df["mean gr"].to_numpy(dtype=float)
da_mean_gr = float(np.mean(da_colony_grs))
da_sem_gr = float(np.std(da_colony_grs, ddof=1) / np.sqrt(len(da_colony_grs)))
da_colony_aucs = da_df["mean auc"].to_numpy(dtype=float)
da_mean_auc = float(np.mean(da_colony_aucs))
da_sem_auc = float(np.std(da_colony_aucs, ddof=1) / np.sqrt(len(da_colony_aucs)))

# Calculate fitness for each colony of each plasmid
# fitness = colony_mean_gr / da28102_mean_gr
order = ["pCU1", "R6K"]
summary = []

for plasmid in order:
    plasmid_df = df.loc[df["plasmid"] == plasmid]
    for _, row in plasmid_df.iterrows():
        colony = int(row["colony"])
        colony_mean_gr = float(row["mean gr"])
        colony_sem_gr = float(row["sem gr"])
        colony_mean_auc = float(row["mean auc"])
        colony_sem_auc = float(row["sem auc"])  
        # fitness = colony growth / da28102 growth
        fitness_gr = colony_mean_gr / da_mean_gr

        # Propagate uncertainty: sigma_f = f * sqrt((sigma_gr/gr)^2 + (sigma_da28102/mg)^2)
        rel_var_gr = (colony_sem_gr / colony_mean_gr) ** 2 if colony_mean_gr > 0 else 0
        rel_var_da = (da_sem_gr / da_mean_gr) ** 2
        fitness_gr_sem = fitness_gr * np.sqrt(rel_var_gr + rel_var_da)

        fitness_auc = colony_mean_auc / da_mean_auc
        
        # Propagate uncertainty: sigma_f = f * sqrt((sigma_gr/gr)^2 + (sigma_da28102    /mg)^2)
        rel_var_auc = (colony_sem_auc / colony_mean_auc) ** 2 if colony_mean_auc > 0 else 0
        rel_var_da = (da_sem_auc / da_mean_auc) ** 2
        fitness_auc_sem = fitness_auc * np.sqrt(rel_var_auc + rel_var_da)
        
        summary.append({
            "plasmid": plasmid,
            "colony": colony,
            "fitness_gr": fitness_gr,
            "fitness_gr_sem": fitness_gr_sem,
            "fitness_auc": fitness_auc,
            "fitness_auc_sem": fitness_auc_sem,
        })

summary = pd.DataFrame(summary)
summary.to_excel("./processed_data/ancestor/conjugative_fitness.xlsx", index=False)

# Summary by plasmid (mean fitness across colonies)
plasmid_summary = []
for plasmid in order:
    plasmid_data = summary.loc[summary["plasmid"] == plasmid]
    mean_fitness_gr = float(np.mean(plasmid_data["fitness_gr"]))
    sem_fitness_gr = float(np.std(plasmid_data["fitness_gr"], ddof=1) / np.sqrt(len(plasmid_data)))
    mean_fitness_auc = float(np.mean(plasmid_data["fitness_auc"]))
    sem_fitness_auc = float(np.std(plasmid_data["fitness_auc"], ddof=1) / np.sqrt(len(plasmid_data)))
    plasmid_summary.append({
        "plasmid": plasmid,
        "mean_fitness_gr": mean_fitness_gr,
        "sem_fitness_gr": sem_fitness_gr,
        "mean_fitness_auc": mean_fitness_auc,
        "sem_fitness_auc": sem_fitness_auc,
    })
plasmid_summary = pd.DataFrame(plasmid_summary)
print(plasmid_summary)
# Plot fitness by plasmid; ax1 = max growth rate, ax2 = AUC
fig, axes = plt.subplots(2,1,figsize=(3,3))
ax1, ax2 = axes
x = np.arange(len(order))
error_kw = {"elinewidth": 1.2, "capsize": 4, "capthick": 1.2}

ax1.bar(x, plasmid_summary["mean_fitness_gr"].to_numpy(dtype=float),
    yerr=plasmid_summary["sem_fitness_gr"].fillna(0).to_numpy(dtype=float),
    width=0.6, error_kw=error_kw, lw=1.2, facecolor="None", edgecolor="k", zorder=2)

# scatters represent fitness for each colony
rng = np.random.default_rng(42)
jitter_scale = 0.08
for j, plasmid in enumerate(order):
    plasmid_data = summary.loc[summary["plasmid"] == plasmid]
    y_vals = plasmid_data["fitness_gr"].to_numpy(dtype=float)
    y_err = plasmid_data["fitness_gr_sem"].to_numpy(dtype=float)
    x_j = rng.normal(loc=x[j], scale=jitter_scale, size=len(y_vals))
    ax1.scatter(x_j, y_vals, s=100, linewidths=1, facecolor=colors[j], edgecolors="k", zorder=3)
    ax1.errorbar(x_j, y_vals, yerr=y_err, marker="o", markersize=0, elinewidth=1.2,
        linewidth=0, capsize=2, capthick=1.2, color="k", zorder=4)
    
ax2.bar(x, plasmid_summary["mean_fitness_auc"].to_numpy(dtype=float),
    yerr=plasmid_summary["sem_fitness_auc"].fillna(0).to_numpy(dtype=float),
    width=0.6, error_kw=error_kw, lw=1.2, facecolor="None", edgecolor="k", zorder=2)

# scatters represent fitness for each colony
rng = np.random.default_rng(42)
jitter_scale = 0.08
for j, plasmid in enumerate(order):
    plasmid_data = summary.loc[summary["plasmid"] == plasmid]
    y_vals = plasmid_data["fitness_auc"].to_numpy(dtype=float)
    y_err = plasmid_data["fitness_auc_sem"].to_numpy(dtype=float)
    x_j = rng.normal(loc=x[j], scale=jitter_scale, size=len(y_vals))
    ax2.scatter(x_j, y_vals, s=100, linewidths=1, facecolor=colors[j], edgecolors="k", zorder=3)
    ax2.errorbar(x_j, y_vals, yerr=y_err, marker="o", markersize=0, elinewidth=1.2,
        linewidth=0, capsize=2, capthick=1.2, color="k", zorder=4)
ax1.axhline(y=1.0, color="k", linestyle="--", linewidth=1, zorder=1)
ax1.set_xticks(x)
ax1.set_xticklabels([])
ax1.set_ylim([0,1.2])
ax1.set_ylabel("w$_{gr}$")
ax2.axhline(y=1.0, color="k", linestyle="--", linewidth=1, zorder=1)
ax2.set_xticks(x)
ax2.set_xticklabels(order, rotation=0)
ax2.set_ylim([0,1.2])
ax2.set_ylabel("w$_{AUC}$")

fig.tight_layout()
fig.savefig("./figures/conjugative_ancestor_fitness.png", dpi=300)
fig.savefig("./figures/conjugative_ancestor_fitness.svg", dpi=300)
plt.show()