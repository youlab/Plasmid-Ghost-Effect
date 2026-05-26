import csv
import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

order = ["pSC101", "colE1", "pUC"]
files = {
    "pSC101": "./processed_data/pSC101_CN.txt",
    "colE1": "./processed_data/colE1_CN.txt",
    "pUC": "./processed_data/pUC_CN.txt",
}

br_means_by_plasmid: dict[str, np.ndarray] = {}
br_sems_by_plasmid: dict[str, np.ndarray] = {}
group_means = []
group_sems = []

for plasmid in order:
    cn = np.loadtxt(files[plasmid])
    mean = np.nanmean(cn, axis=1)
    std = np.nanstd(cn, axis=1, ddof=1)
    sem = std / np.sqrt(cn.shape[1])

    br_means_by_plasmid[plasmid] = mean
    br_sems_by_plasmid[plasmid] = sem

    group_means.append(float(np.nanmean(mean)))
    group_sems.append(float(np.nanstd(mean, ddof=1) / np.sqrt(len(mean))))

group_means = np.array(group_means, dtype=float)
group_sems = np.array(group_sems, dtype=float)

report_path = "./processed_data/PCN_summary_report.tsv"
with open(report_path, "w", newline="") as report_file:
    writer = csv.writer(report_file, delimiter="\t")
    writer.writerow(["plasmid", "mean", "sem"])
    for plasmid, mean_value, sem_value in zip(order, group_means, group_sems):
        writer.writerow([plasmid, f"{mean_value:.6f}", f"{sem_value:.6f}"])

fig_pcn, ax_pcn = plt.subplots(figsize=(3, 3))
x = np.arange(len(order))

error_kw = {"elinewidth": 1.2, "capsize": 4, "capthick": 1.2}
ax_pcn.bar(
    x,
    group_means,
    yerr=group_sems,
    width=0.4,
    error_kw=error_kw,
    lw=1.2,
    facecolor="None",
    edgecolor="k",
    zorder=2,
)

# Scatters represent per-colony mean ± sem from technical replicates.
rng = np.random.default_rng(42)
jitter_scale = 0.08
for j, plasmid in enumerate(order):
    y_vals = br_means_by_plasmid[plasmid]
    y_err = br_sems_by_plasmid[plasmid]
    x_j = rng.normal(loc=x[j], scale=jitter_scale, size=len(y_vals))

    ax_pcn.scatter(
        x_j,
        y_vals,
        s=60,
        linewidths=1,
        facecolor=colors[j],
        edgecolors="k",
        zorder=3,
    )
    ax_pcn.errorbar(
        x_j,
        y_vals,
        yerr=y_err,
        marker="o",
        markersize=0,
        elinewidth=1.2,
        linewidth=0,
        capsize=2,
        capthick=1.2,
        color="k",
        zorder=4,
    )

ax_pcn.set_xticks(x)
ax_pcn.set_xticklabels(order, rotation = 90)
ax_pcn.set_ylabel("plasmid copy number")
fig_pcn.tight_layout()

fig_pcn.savefig("./figures/PCN_summary.png", dpi=300)
fig_pcn.savefig("./figures/PCN_summary.svg", dpi=300)
plt.show()