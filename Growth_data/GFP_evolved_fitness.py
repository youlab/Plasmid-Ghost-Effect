import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from calculate_growth_rate import GR

day = 3
plasmid = "pUC"

cmap = plt.get_cmap('Set2')
colors = ["#AAAAAA"]+[cmap(i) for i in np.linspace(0, 1, 8)]
plasmid_id = {"pSC101":1, "colE1":2, "pUC":3}[plasmid]
plasmid_color = colors[plasmid_id]

processed_dir = "./processed_data"

def load_od_curve(xlsx_path: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return time, replicate matrix, and replicate labels."""
    df = pd.read_excel(xlsx_path)
    time = df["time"].to_numpy(dtype=float)
    replicate_cols = [c for c in df.columns if c != "time"]
    od_values = df[replicate_cols].to_numpy(dtype=float)
    return time, od_values, replicate_cols

# AUC calculation from saved OD xlsx (one value per replicate)
def quantify_fitness(time: np.ndarray, od_matrix: np.ndarray, labels: list[str], plasmid: str, biorep: int) -> pd.DataFrame:
    records = []
    for colony in range(1,4,1):
        auc = []
        gr = []
        rep_labels = [label for label in labels if re.fullmatch(rf"c{colony}_(\d+)", label)]
        for rep_label in rep_labels:
            od = od_matrix[:, labels.index(rep_label)]
            auc.append(float(np.sum((od[:-1] + od[1:]) * np.diff(time) * 0.5)))
            od = np.maximum(od, 0.01)  # floor to avoid log(0) and small value fluctuations
            gr.append(np.max(GR(time[7:], np.log(od[7:]), window_size=6)[1]))  # max growth rate from rolling window
        if not rep_labels:
            continue
        records.append(
            {
                "plasmid": plasmid,
                "time": day,
                "replicate": biorep,
                "colony": colony,
                "mean auc": np.mean(auc),
                "sem auc": np.std(auc, ddof=1) / np.sqrt(len(auc)),
                "mean gr": np.mean(gr),
                "sem gr": np.std(gr, ddof=1) / np.sqrt(len(gr))
            }
        )
    return pd.DataFrame(records)

fig, axes = plt.subplots(3, 3, figsize=(6, 6))
df = []
for br in range(3):
    time_mg, mg_od, mg_reps = load_od_curve(f"./processed_data/Day{day}/{plasmid}-BioRep{br+1}_OD.xlsx")
    time_plasmid, plasmid_od, plasmid_reps = load_od_curve(f"./processed_data/Day{day}/{plasmid}+BioRep{br+1}_OD.xlsx")
    df.append(quantify_fitness(time_mg, mg_od, mg_reps, f"{plasmid}-", br+1))
    df.append(quantify_fitness(time_plasmid, plasmid_od, plasmid_reps, f"{plasmid}+", br+1))
    for col in range(3):
        ax = axes[br, col]
        rep_label = [label for label in mg_reps if re.fullmatch(rf"c{col+1}_(\d+)", label)]
        mg_idx = [mg_reps.index(label) for label in rep_label]
        od = mg_od[:, mg_idx]
        od = np.maximum(od, 0.01)
        mg_mean = np.mean(od, axis=1)
        mg_std = np.std(od, axis=1, ddof=1)
        t, gr = GR(time_mg[7:], np.log(mg_mean[7:]), window_size=6)
        t_grmax = t[np.argmax(gr)]
        y_grmax = mg_mean[int(t_grmax * 12)]
        ax.plot(time_mg, mg_mean, linewidth=1.5, c=colors[0], label=f"{plasmid}-")
        ax.fill_between(time_mg, mg_mean - mg_std, mg_mean + mg_std, color=colors[0], alpha=0.25)
        ax.scatter(t_grmax, y_grmax, fc='None', ec=colors[0], lw = 2, s=50, zorder=5)

        rep_label = [label for label in plasmid_reps if re.fullmatch(rf"c{col+1}_(\d+)", label)]
        plasmid_idx = [plasmid_reps.index(label) for label in rep_label]
        od = plasmid_od[:, plasmid_idx]
        od = np.maximum(od, 0.01)
        plasmid_mean = np.mean(od, axis=1)
        plasmid_std = np.std(od, axis=1, ddof=1)
        t, gr = GR(time_plasmid[7:], np.log(plasmid_mean[7:]), window_size=6)
        t_grmax = t[np.argmax(gr)]
        y_grmax = plasmid_mean[int(t_grmax * 12)] 
        ax.plot(time_plasmid, plasmid_mean, linewidth=1.5, c=plasmid_color, label=f"{plasmid}+")
        ax.fill_between(time_plasmid, plasmid_mean - plasmid_std, plasmid_mean + plasmid_std, color=plasmid_color, alpha=0.25)
        ax.scatter(t_grmax, y_grmax, fc='None', ec=plasmid_color, lw = 2, s=50, zorder=5)

        ax.text(x=0.1,y=0.9,s=f"bio rep {br+1}\ncolony {col+1}", fontsize=12, va="top", ha="left", transform=ax.transAxes)
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_xticks([0, 12, 24])
        ax.set_ylim([0,1])
        
        if br == 2 and col == 0:
            ax.set_xlabel("time (hours)")
            ax.set_ylabel("OD600")
        # else:
        #     ax.set_xticklabels([])
        #     ax.set_yticklabels([])
        if br == 0 and col == 1:
            ax.set_title(f"MG1655 - {plasmid} Day {day}", fontsize=14)
        if br == 0 and col == 0:
            ax.legend(
                loc="lower right",
                fontsize=12,
                handlelength=0.8,
                handletextpad=0.4,
                labelspacing=0.35,
                frameon=False,
            )
fig.tight_layout()
fig.savefig(f"./figures/{plasmid}_Day{day}_growth_curve.png", dpi=300)

df = pd.concat(df, ignore_index=True)

# Calculate fitness for each colony of each biological replicate
# fitness = plasmid+ colony_gr / plasmid- colony_gr (and auc as well)
summary = []

for br in range(1, 4):
    for colony in range(1, 4):
        # Get plasmid- (cured) data
        cured_row = df.loc[(df["replicate"] == br) & (df["colony"] == colony) & (df["plasmid"] == f"{plasmid}-")]
        # Get plasmid+ (carrying) data
        carrying_row = df.loc[(df["replicate"] == br) & (df["colony"] == colony) & (df["plasmid"] == f"{plasmid}+")]
        
        if len(cured_row) == 0 or len(carrying_row) == 0:
            continue
        
        cured_gr = float(cured_row["mean gr"].iloc[0])
        cured_gr_sem = float(cured_row["sem gr"].iloc[0])
        cured_auc = float(cured_row["mean auc"].iloc[0])
        cured_auc_sem = float(cured_row["sem auc"].iloc[0])
        
        carrying_gr = float(carrying_row["mean gr"].iloc[0])
        carrying_gr_sem = float(carrying_row["sem gr"].iloc[0])
        carrying_auc = float(carrying_row["mean auc"].iloc[0])
        carrying_auc_sem = float(carrying_row["sem auc"].iloc[0])
        
        # Fitness = carrying / cured
        fitness_gr = carrying_gr / cured_gr
        # Propagate uncertainty: sigma_f = f * sqrt((sigma_carrying/carrying)^2 + (sigma_cured/cured)^2)
        rel_var_gr = (carrying_gr_sem / carrying_gr) ** 2
        rel_var_cured_gr = (cured_gr_sem / cured_gr) ** 2
        fitness_gr_sem = fitness_gr * np.sqrt(rel_var_gr + rel_var_cured_gr)
        
        fitness_auc = carrying_auc / cured_auc
        # Propagate uncertainty
        rel_var_auc = (carrying_auc_sem / carrying_auc) ** 2
        rel_var_cured_auc = (cured_auc_sem / cured_auc) ** 2
        fitness_auc_sem = fitness_auc * np.sqrt(rel_var_auc + rel_var_cured_auc)
        
        summary.append({
            "day": day,
            "plasmid": plasmid,
            "replicate": br,
            "colony": colony,
            "fitness_gr": fitness_gr,
            "fitness_gr_sem": fitness_gr_sem,
            "fitness_auc": fitness_auc,
            "fitness_auc_sem": fitness_auc_sem,
        })

summary = pd.DataFrame(summary)
summary.to_excel(f"./processed_data/Day{day}/{plasmid}_fitness.xlsx", index=False)
plt.show()