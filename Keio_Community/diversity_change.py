import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def p_to_star(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

communities = ["Comm87","Comm57_R6K_100","Comm57_pCU1_100","Comm57_R388_100","Comm57_RP4_100",
               "Comm57_R6K_500","Comm57_pCU1_500","Comm57_R388_500","Comm57_RP4_500"]
plasmids = ["R6K","pCU1","R388","RP4"]
colors = ['#8da0cb', '#fc8d62']

Labels=["no pulse","+Ab",]

fig, axes = plt.subplots(2, 5, figsize=(10, 4))
axes=axes.flat
for i, comm in enumerate(communities):
    data = np.load(f"./composition_all/{comm}.npy")
    time = np.array([0,2,3,7,10,15,20,25,30,35]) if i==0 else np.array([0,2,3,7,10,15,20])
    simpson = np.zeros((6, len(time)))
    if i <5:
        ax=axes[i]
    else:
        ax = axes[i + 1]

    for j in range(6):
        comp = data[j,:]
        div = 1/np.sum(comp**2,axis=0)
        simpson[j,:]=div
        k = j//3
        ax.plot(time, div, color=colors[k], linewidth=1.5,zorder=0)
        ax.scatter(time, div, s=30, color=colors[k], linewidth=0.5, edgecolors='black',zorder=10)
    ax.fill_between(x=[2, 3], y1=[0.5, 0.5], y2=[80,80],
                    color="#808080", zorder=-20, alpha=0.3)
    ax.set_yscale("log")
    div_np = simpson[0:3,2]
    div_p = simpson[3:,2]
    t_stat, p_value = stats.ttest_ind(div_np, div_p, equal_var=False) # two-sided Welch's t-test
    print(f"{comm} Welch's t-test (two-sided): {p_value:.2e}")
    ax.text(3, 50, p_to_star(p_value), ha='center', va='top', fontsize=12)
    if i==0:
        ax.set_xticks([0,15,30])
        ax.text(x=0.9,y=0.5,s="Comm87\nR388",ha="right",transform=ax.transAxes)
    else:
        ax.set_xticks([0, 10, 20])
        if i<5:
            plasmid = plasmids[i%5-1]
            ratio = 100
        else:
            plasmid = plasmids[i % 5]
            ratio = 500
        if i in [1,5]:
            ax.text(x=0.9, y=0.2, s=f"Comm57\n{plasmid}\n1:{ratio}", ha="right", transform=ax.transAxes)
        else:
            ax.text(x=0.9, y=0.05, s=f"{plasmid}\n1:{ratio}", ha="right", transform=ax.transAxes)
    ax.set_yticks([1,10])
    ax.set_ylim([0.5,60])
axes[5].scatter([], [], s=30, color=colors[0], linewidth=0.5, edgecolors='black',label="no pulse")
axes[5].scatter([], [], s=30, color=colors[1], linewidth=0.5, edgecolors='black',label="+ Ab")
axes[5].set_visible(False)
fig.legend(loc="lower left",bbox_to_anchor=(0.08,0.2))
axes[0].set_xlabel("time (days)")
axes[0].set_ylabel("$^2$D",rotation=0,va="center",ha="right")
fig.subplots_adjust(bottom=0.1,left=0.08,top=0.95,right=0.98,hspace=0.25,wspace=0.25)
fig.savefig("./figures/simpson_diversity.png",dpi=300)

plt.show()
