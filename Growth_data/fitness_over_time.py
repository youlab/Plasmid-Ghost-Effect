import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy import stats

cmap = plt.get_cmap('Set2')
colors = ["#AAAAAA"]+[cmap(i) for i in np.linspace(0, 1, 8)]
plasmid_id = {"pSC101":1, "colE1":2, "pUC":3}

# Containers to collect p-values computed in the plotting loop
collected_gr_pvals = {}
collected_auc_pvals = {}

def add_significance_brackets(ax, all_y_vals, pvals_dict):
    """Add significance brackets to bar plot with proper layering.
    
    Parameters:
    - all_y_vals: list of arrays, where all_y_vals[i] contains all y values for source i
    - pvals_dict: dict mapping (i_a, i_b) to p-value
    """
    # Find the maximum y value across all data points
    max_vals = [np.nanmax(arr) if len(arr) > 0 else 0 for arr in all_y_vals]
    y_max = np.nanmax(max_vals) if max_vals else 1
    bracket_spacing = 0.13
    y_offset = 0.1
    # Sort by span distance to layer them properly
    sorted_comparisons = sorted(pvals_dict.items(), key=lambda x: x[0][1] - x[0][0])
    
    layer_map = {}
    layer = -1
    for idx, ((i_a, i_b), p_val) in enumerate(sorted_comparisons):
        if p_val < 0.05:
            layer += 1
        layer_map[idx] = layer
    
    for idx, ((i_a, i_b), p_val) in enumerate(sorted_comparisons):
        if p_val < 0.05:
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*"
            layer = layer_map[idx]
            y_bracket = y_max + y_offset + layer * bracket_spacing
            ax.plot([i_a, i_b], [y_bracket,y_bracket], 'k-', linewidth=1, zorder=5)
            # Add significance label
            ax.text((i_a + i_b) / 2, y_bracket, sig, 
                   ha='center', fontsize=14, zorder=5)

# Load ancestor fitness
ancestor_fitness = pd.read_excel("./processed_data/ancestor/GFP_fitness.xlsx")
ancestor_fitness["source"] = "ancestor"
ancestor_fitness["day"] = np.nan

# Load evolved fitness for all days and plasmids
plasmids = ["pSC101", "colE1", "pUC"]
days = [3, 20]

evolved_fitness_list = []
for day in days:
    for plasmid in plasmids:
        filepath = f"./processed_data/Day{day}/{plasmid}_fitness.xlsx"
        if os.path.exists(filepath):
            df = pd.read_excel(filepath)
            df["plasmid"] = plasmid
            df["day"] = day
            df["source"] = f"Day{day}"
            evolved_fitness_list.append(df)
        else:
            print(f"Warning: File not found: {filepath}")

evolved_fitness = pd.concat(evolved_fitness_list, ignore_index=True)

# Create combined dataframe with consistent columns
combined = pd.concat([ancestor_fitness, evolved_fitness], ignore_index=True)
combined.to_excel("./processed_data/combined_fitness.xlsx", index=False)
print("Saved combined fitness to ./processed_data/combined_fitness.xlsx")

# Create visualizations by plasmid
fig_main, axes_main = plt.subplots(3, 2, figsize=(8, 10))

for i, plasmid in enumerate(plasmids):
    plasmid_data = combined.loc[combined["plasmid"] == plasmid]
    
    # GR plot
    ax_gr = axes_main[i, 0]
    # AUC plot
    ax_auc = axes_main[i, 1]
    
    # Prepare data for plotting
    sources = ["ancestor"] + [f"Day{d}" for d in days]
    gr_means = []
    gr_sems = []
    auc_means = []
    auc_sems = []
    
    for source in sources:
        source_data = plasmid_data.loc[plasmid_data["source"] == source]
        if len(source_data) > 0:
            gr_means.append(float(np.mean(source_data["fitness_gr"])))
            gr_sems.append(float(np.std(source_data["fitness_gr"], ddof=1) / np.sqrt(len(source_data))))
            auc_means.append(float(np.mean(source_data["fitness_auc"])))
            auc_sems.append(float(np.std(source_data["fitness_auc"], ddof=1) / np.sqrt(len(source_data))))
        else:
            gr_means.append(np.nan)
            gr_sems.append(np.nan)
            auc_means.append(np.nan)
            auc_sems.append(np.nan)
    
    x = np.arange(len(sources))
    error_kw = {"elinewidth": 1.2, "capsize": 4, "capthick": 1.2}
    
    # Pre-compute t-tests for GR and AUC
    gr_pvals = {}
    auc_pvals = {}
    for i_a, source_a in enumerate(sources):
        for i_b, source_b in enumerate(sources[i_a+1:], start=i_a+1):
            data_a = plasmid_data.loc[plasmid_data["source"] == source_a, "fitness_gr"].dropna().to_numpy()
            data_b = plasmid_data.loc[plasmid_data["source"] == source_b, "fitness_gr"].dropna().to_numpy()
            if len(data_a) > 0 and len(data_b) > 0:
                _, p = stats.ttest_ind(data_a, data_b, equal_var=False)
                gr_pvals[(i_a, i_b)] = p
            
            data_a = plasmid_data.loc[plasmid_data["source"] == source_a, "fitness_auc"].dropna().to_numpy()
            data_b = plasmid_data.loc[plasmid_data["source"] == source_b, "fitness_auc"].dropna().to_numpy()
            if len(data_a) > 0 and len(data_b) > 0:
                _, p = stats.ttest_ind(data_a, data_b, equal_var=False)
                auc_pvals[(i_a, i_b)] = p
    
    # GR bar plot
    ax_gr.bar(x, gr_means, yerr=gr_sems, width=0.6, error_kw=error_kw, 
              lw=1.2, facecolor="None", edgecolor="k", zorder=2)
    
    # GR scatter for individual colonies/replicates
    rng = np.random.default_rng(42)
    jitter_scale = 0.08
    gr_all_y_vals = []
    for j, source in enumerate(sources):
        source_data = plasmid_data.loc[plasmid_data["source"] == source]
        if len(source_data) > 0:
            y_vals = source_data["fitness_gr"].to_numpy(dtype=float)
            gr_all_y_vals.append(y_vals)
            y_err = source_data.get("fitness_gr_sem", pd.Series([0]*len(y_vals))).to_numpy(dtype=float)
            x_j = rng.normal(loc=x[j], scale=jitter_scale, size=len(y_vals))
            ax_gr.scatter(x_j, y_vals, s=60, linewidths=1, facecolor=colors[plasmid_id[plasmid]], 
                         edgecolors="k", zorder=3)
            if len(y_err) > 0 and np.any(y_err > 0):
                ax_gr.errorbar(x_j, y_vals, yerr=y_err, marker="o", markersize=0, elinewidth=1.2,
                    linewidth=0, capsize=2, capthick=1.2, color="k", zorder=4)
        else:
            gr_all_y_vals.append(np.array([]))
    
    # Add significance brackets for GR
    add_significance_brackets(ax_gr, gr_all_y_vals, gr_pvals)
    
    # AUC bar plot
    ax_auc.bar(x, auc_means, yerr=auc_sems, width=0.6, error_kw=error_kw,
               lw=1.2, facecolor="None", edgecolor="k", zorder=2)
    
    # AUC scatter for individual colonies/replicates
    rng = np.random.default_rng(42)
    auc_all_y_vals = []
    for j, source in enumerate(sources):
        source_data = plasmid_data.loc[plasmid_data["source"] == source]
        if len(source_data) > 0:
            y_vals = source_data["fitness_auc"].to_numpy(dtype=float)
            auc_all_y_vals.append(y_vals)
            y_err = source_data.get("fitness_auc_sem", pd.Series([0]*len(y_vals))).to_numpy(dtype=float)
            x_j = rng.normal(loc=x[j], scale=jitter_scale, size=len(y_vals))
            ax_auc.scatter(x_j, y_vals, s=60, linewidths=1, facecolor=colors[plasmid_id[plasmid]],
                          edgecolors="k", zorder=3)
            if len(y_err) > 0 and np.any(y_err > 0):
                ax_auc.errorbar(x_j, y_vals, yerr=y_err, marker="o", markersize=0, elinewidth=1.2,
                    linewidth=0, capsize=2, capthick=1.2, color="k", zorder=4)
        else:
            auc_all_y_vals.append(np.array([]))
    
    # Add significance brackets for AUC
    add_significance_brackets(ax_auc, auc_all_y_vals, auc_pvals)
    # Save p-values for later reporting (do not re-run tests later)
    collected_gr_pvals[plasmid] = gr_pvals
    collected_auc_pvals[plasmid] = auc_pvals
    
    # Format GR subplot
    ax_gr.axhline(y=1.0, color="k", linestyle="--", linewidth=1, zorder=1)
    ax_gr.set_xticks(x)
    ax_gr.set_xticklabels(sources, rotation=45, ha="right")
    ax_gr.set_ylim([0, 1.6])
    ax_gr.set_yticks([0, 0.5, 1.0, 1.5])
    ax_gr.set_ylabel("w$_{gr}$", fontsize=12)
    ax_gr.set_title(f"{plasmid} - Growth Rate", fontsize=12)
    
    # Format AUC subplot
    ax_auc.axhline(y=1.0, color="k", linestyle="--", linewidth=1, zorder=1)
    ax_auc.set_xticks(x)
    ax_auc.set_xticklabels(sources, rotation=45, ha="right")
    ax_auc.set_ylim([0, 1.6])
    ax_auc.set_yticks([0, 0.5, 1.0, 1.5])
    ax_auc.set_ylabel("w$_{AUC}$", fontsize=12)
    ax_auc.set_title(f"{plasmid} - AUC", fontsize=12)

fig_main.tight_layout()
fig_main.savefig("./figures/combined_fitness_all_plasmids.png", dpi=300)
#fig_main.savefig("./figures/combined_fitness_all_plasmids.svg", dpi=300)
print("Saved combined fitness plots")

# Create summary table by plasmid and day
summary_rows = []
for plasmid in plasmids:
    plasmid_data = combined.loc[combined["plasmid"] == plasmid]
    for source in ["ancestor"] + [f"Day{d}" for d in days]:
        source_data = plasmid_data.loc[plasmid_data["source"] == source]
        if len(source_data) > 0:
            summary_rows.append({
                "plasmid": plasmid,
                "timepoint": source,
                "mean_fitness_gr": float(np.mean(source_data["fitness_gr"])),
                "sem_fitness_gr": float(np.std(source_data["fitness_gr"], ddof=1) / np.sqrt(len(source_data))),
                "mean_fitness_auc": float(np.mean(source_data["fitness_auc"])),
                "sem_fitness_auc": float(np.std(source_data["fitness_auc"], ddof=1) / np.sqrt(len(source_data))),
                "n_samples": len(source_data),
            })

summary = pd.DataFrame(summary_rows)
#summary.to_excel("./processed_data/fitness_summary_by_plasmid_timepoint.xlsx", index=False)
#print("Saved fitness summary table")
print("\nFitness Summary:")
print(summary.to_string(float_format="{:.2f}".format))

# Print stored p-values computed earlier during plotting (no re-running of stats)
print("\n" + "="*80)
print("PRECOMPUTED PAIRWISE P-VALUES (from plotting step)")
print("="*80)
source_names = ["ancestor", "Day3", "Day20"]
for plasmid in plasmids:
    print(f"\n{plasmid}:")
    # GR
    print("  Growth Rate (fitness_gr):")
    gr_p = collected_gr_pvals.get(plasmid, {})
    if not gr_p:
        print("    (no comparisons computed)")
    for (i_a, i_b), p in sorted(gr_p.items()):
        s1, s2 = source_names[i_a], source_names[i_b]
        d1 = combined.loc[(combined["plasmid"] == plasmid) & (combined["source"] == s1), "fitness_gr"].dropna().to_numpy()
        d2 = combined.loc[(combined["plasmid"] == plasmid) & (combined["source"] == s2), "fitness_gr"].dropna().to_numpy()
        n1, n2 = len(d1), len(d2)
        mean1 = np.mean(d1) if n1>0 else np.nan
        mean2 = np.mean(d2) if n2>0 else np.nan
        sd1 = np.std(d1, ddof=1) if n1>1 else np.nan
        sd2 = np.std(d2, ddof=1) if n2>1 else np.nan
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"    {s1} (n={n1}, μ={mean1:.2f}±{sd1:.2f}) vs {s2} (n={n2}, μ={mean2:.2f}±{sd2:.2f}): p={p:.3f} {sig}")
    # AUC
    print("  AUC (fitness_auc):")
    auc_p = collected_auc_pvals.get(plasmid, {})
    if not auc_p:
        print("    (no comparisons computed)")
    for (i_a, i_b), p in sorted(auc_p.items()):
        s1, s2 = source_names[i_a], source_names[i_b]
        d1 = combined.loc[(combined["plasmid"] == plasmid) & (combined["source"] == s1), "fitness_auc"].dropna().to_numpy()
        d2 = combined.loc[(combined["plasmid"] == plasmid) & (combined["source"] == s2), "fitness_auc"].dropna().to_numpy()
        n1, n2 = len(d1), len(d2)
        mean1 = np.mean(d1) if n1>0 else np.nan
        mean2 = np.mean(d2) if n2>0 else np.nan
        sd1 = np.std(d1, ddof=1) if n1>1 else np.nan
        sd2 = np.std(d2, ddof=1) if n2>1 else np.nan
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"    {s1} (n={n1}, μ={mean1:.2f}±{sd1:.2f}) vs {s2} (n={n2}, μ={mean2:.2f}±{sd2:.2f}): p={p:.1e} {sig}")

plt.show()
