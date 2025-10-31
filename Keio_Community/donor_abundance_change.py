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
p_values = []
delta_max = []
donor_ids = [0,3,2,0,1]

for i, comm in enumerate(communities):
    donor = donor_ids[i]
    data = np.load(f"./composition_all/{comm}.npy")
    donor_day2 = data[:,donor,1]
    donor_day3 = data[:,donor,2]
    delta = donor_day3-donor_day2
    delta_NoAb = delta[0:3]
    delta_Ab = delta[3:]
    delta_max.append(np.max(np.abs(delta)))
    t_stat, p_value = stats.ttest_ind(delta_NoAb, delta_Ab, equal_var=False) # two-sided Welch's t-test
    print(f"{comm} Welch's t-test (two-sided): {p_value:.2e}")
    p_values.append(p_value)
    df["delta"]+=list(delta)
    df["condition"]+=(["no pulse"]*3+["+Ab"]*3)
    df["community"]+=[comm]*6
df=pd.DataFrame(df)
fig, ax = plt.subplots(1, 1, figsize=(3, 1.9))
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
fig.savefig("./figures/donor_abundance_change.png",dpi=300)
fig.savefig("./figures/donor_abundance_change.svg")
plt.show()