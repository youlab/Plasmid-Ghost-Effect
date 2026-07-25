import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# --- statistics on the changes in donor abundance ---
colors = ['#8da0cb', '#fc8d62']

def p_to_star(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

df = {"delta":[],"condition":[],"community":[]}
communities = ["Comm87","Comm57_R6K_100","Comm57_pCU1_100","Comm57_R388_100","Comm57_RP4_100"]
plasmids = ["R388","R6K","pCU1","R388","RP4"]
p_values = []
delta_max = []
donor_ids = [0,3,2,0,1]
stats_rows = []

for i, comm in enumerate(communities):
    donor = donor_ids[i]
    data = np.load(f"./composition_all/{comm}.npy")
    donor_day2 = data[:,donor,1]
    donor_day3 = data[:,donor,2]
    delta = donor_day3-donor_day2
    delta_NoAb = delta[0:3]
    delta_Ab = delta[3:]
    delta_max.append(np.max(np.abs(delta)))
    # contrast is +Ab - no pulse, so the sign matches the direction of the antibiotic effect
    t_stat, p_value = stats.ttest_ind(delta_Ab, delta_NoAb, equal_var=False) # two-sided Welch's t-test
    n1, n2 = len(delta_NoAb), len(delta_Ab)
    s1, s2 = delta_NoAb.std(ddof=1), delta_Ab.std(ddof=1)
    # Welch-Satterthwaite degrees of freedom
    v_NoAb = s1**2/n1
    v_Ab = s2**2/n2
    dof = (v_NoAb+v_Ab)**2/(v_NoAb**2/(n1-1)+v_Ab**2/(n2-1))
    # 95% CI of the difference of means (Welch)
    mean_diff = delta_Ab.mean()-delta_NoAb.mean()
    se_diff = np.sqrt(v_NoAb+v_Ab)
    t_crit = stats.t.ppf(0.975, dof)
    ci_low, ci_high = mean_diff-t_crit*se_diff, mean_diff+t_crit*se_diff
    # effect size: Cohen's d (pooled sd), 95% CI from the normal-approximation SE
    s_pooled = np.sqrt(((n1-1)*s1**2+(n2-1)*s2**2)/(n1+n2-2))
    d = mean_diff/s_pooled
    se_d = np.sqrt((n1+n2)/(n1*n2)+d**2/(2*(n1+n2-2)))
    d_low, d_high = d-1.96*se_d, d+1.96*se_d
    p_values.append(p_value)
    stats_rows.append({
        "comparison": f"{plasmids[i]} in {comm.split('_')[0]}",
        "n_per_group": n1,
        "mean_no_pulse": delta_NoAb.mean(), "sd_no_pulse": s1,
        "mean_Ab": delta_Ab.mean(), "sd_Ab": s2,
        "mean_diff": mean_diff, "diff_CI_low": ci_low, "diff_CI_high": ci_high,
        "cohens_d": d, "d_CI_low": d_low, "d_CI_high": d_high,
        "t": t_stat, "df": dof, "p": p_value, "sig": p_to_star(p_value),
    })
    df["delta"]+=list(delta)
    df["condition"]+=(["no pulse"]*3+["+Ab"]*3)
    df["community"]+=[comm]*6
df=pd.DataFrame(df)
df.to_csv("./figures/donor_abundance_change_data.csv", index=False) # per-replicate source data

stats_table = pd.DataFrame(stats_rows)
stats_table.to_csv("./figures/donor_abundance_change_stats.csv", index=False) # full precision

# printed table: n, t, p, mean difference [95% CI], d -- rounded to one decimal (p exact)
report = pd.DataFrame({
    "comparison": stats_table["comparison"],
    "n": stats_table["n_per_group"],
    "t (df)": [f"{r.t:.1f} ({r.df:.1f})" for r in stats_table.itertuples()],
    "p": [f"{r.p:.3e} {r.sig}" for r in stats_table.itertuples()],
    "mean diff [95% CI]": [f"{r.mean_diff:.2f} [{r.diff_CI_low:.2f}, {r.diff_CI_high:.2f}]"
                           for r in stats_table.itertuples()],
    "Cohen's d [95% CI]": [f"{r.cohens_d:.1f} [{r.d_CI_low:.1f}, {r.d_CI_high:.1f}]"
                           for r in stats_table.itertuples()],
})
print("--- Delta donor: +Ab vs no pulse, two-sided Welch's t-test ---")
print(report.to_string(index=False))
fig, ax = plt.subplots(1, 1, figsize=(4, 1.9))
error_kw = dict(linewidth=0.8, color="k")
sns.barplot(data=df,x="community",y="delta",hue="condition",ax=ax,errorbar="sd",facecolor="None",edgecolor="k",err_kws=error_kw, capsize=0.2, lw=0.8, legend=False)
sns.stripplot(data=df,x="community",y="delta",hue="condition",ax=ax,palette=colors,dodge=True,jitter=True,zorder=-10)
for i,p in enumerate(p_values):
    ymax = delta_max[i]
    if i!=4:
        ax.text(i, ymax + 0.05, p_to_star(p), ha='center', va='bottom', fontsize=12)
        ax.plot([i-0.2,i+0.2], np.ones(2) * ymax + 0.2, c="k", lw=1)
    else:
        ax.text(i, - ymax - 0.25, p_to_star(p), ha='center', va='top', fontsize=12)
        ax.plot([i-0.2,i+0.2], - np.ones(2) * ymax - 0.2, c="k", lw=1)
ax.legend(loc="lower left",bbox_to_anchor=(0,-0.05),labelspacing=0.1,frameon=False,handlelength=0.2)
ax.set_xticks(np.arange(0,5,1))
ax.set_xticklabels(["R388","R6K","pCU1","R388","RP4"],rotation=20)
ax.set_xlabel("")
ax.set_ylabel(r"$\Delta$ donor")
ax.set_yticks([-1,0,1])
ax.set_ylim([-1.5,1.5])
xlim = ax.get_xlim()
ax.plot(xlim,np.zeros(2),c="k",lw=0.8)
ax.set_xlim(xlim)
fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.87)
# fig.savefig("./figures/donor_abundance_change.png",dpi=300)
# fig.savefig("./figures/donor_abundance_change.svg")
plt.show()