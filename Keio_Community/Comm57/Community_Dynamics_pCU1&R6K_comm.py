import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plasmid = "R6K"
Ab = "Carb" if plasmid == "pCU1" else "Strp"
dilution = 100

time=np.array([0,2,3,7,10,15,20])
Conditions=["NoAb","Ab",]
Labels=["no pulse",f"+{Ab}",]
Reps=3
n_show = 12
cmap = plt.get_cmap('Set3')
colors = [cmap(i) for i in np.linspace(0, 1, n_show)] + ["#B3B3B3"]

all_data = []
for i in range(len(Conditions)):
    for j in range(Reps):
        df = pd.read_excel(f"./processed_data/organized_NGS/{plasmid}_{dilution}_dilution.xlsx",
                            sheet_name=f"{Conditions[i]}_Rep{j+1}",index_col=0)
        df=np.array(df)[:-2]#exclude rows "total" and "unmatched
        comp = df/np.sum(df,axis=0)
        all_data.append(comp)
all_data = np.array(all_data)

sum = np.sum(all_data,axis=(0,2))
rank_idx = np.argsort(-sum)
donors = np.arange(4)
rank_idx = np.r_[donors, rank_idx[~np.isin(rank_idx, donors)]] # move the donor strains (1-4) to the front
strain_ids = pd.read_excel(f"./processed_data/organized_NGS/{plasmid}_{dilution}_dilution.xlsx",
                           sheet_name=f"NoAb_Rep1")["Barcodes"]
strain_ids = np.array(strain_ids)[:-2]
strain_ids = strain_ids[rank_idx]
print(strain_ids)

for i in range(len(Conditions)):
    fig, ax = plt.subplots(1, 1, figsize=(2.05, 1.9))
    composition = []
    for j in range(Reps):
        df = pd.read_excel(f"./processed_data/organized_NGS/{plasmid}_{dilution}_dilution.xlsx",
                            sheet_name=f"{Conditions[i]}_Rep{j+1}",index_col=0)
        df = np.array(df)[:-2]  # exclude rows "total" and "unmatched
        comp = df/np.sum(df,axis=0)
        composition.append(comp)
    composition = np.array(composition)
    composition = composition[:,rank_idx,:]
    mean_composition = np.mean(composition,axis=0)
    std_composition = np.std(composition,axis=0,ddof=1)
    n = mean_composition.shape[0]
    base = np.zeros(len(time))
    for k in range(n_show):
        mean = mean_composition[k,:]
        std = std_composition[k,:]
        se = std/np.sqrt(3)
        ax.errorbar(time, base+mean, se, marker='None', linewidth=0.5,
                    elinewidth=0.5, capsize=2, capthick=0.5, color="k")
        ax.fill_between(time, base, base+mean, color=colors[k])
        base = base+mean
    top_comp = np.sum(mean_composition[0:n_show,:],axis=0)
    ax.fill_between(time, top_comp, np.ones_like(time), color=colors[-1])
    ax.set_yticks([])
    ax.set_ylim([0,1])
    ax.set_xticks(np.arange(0,21,5))
    ax.set_xlim([0,20])
    ax.set_title(Labels[i])

    if i==0:
        #ax.set_ylabel("Comm57")
        ax.set_xlabel("time (days)")
    fig.subplots_adjust(left=0.25, right=0.95, bottom=0.15, top=0.87)
    fig.savefig(f"./figures/LT21_comm_dynamics_{plasmid}_{Conditions[i]}.png",dpi=300)
    fig.savefig(f"./figures/LT21_comm_dynamics_{plasmid}_{Conditions[i]}.svg")

from matplotlib.patches import Patch
fig2,ax2=plt.subplots(1,1)
labels = list(strain_ids[0:n_show])+["Other"]
legend_handles = [
    Patch(facecolor=colors[i], label=labels[i])
    for i in range(2, n_show+1) # exclude strain 1 and 2 since they are not present
]

ax2.legend(handles=legend_handles,title="strain id")
fig2.savefig(f"./figures/LT21_legends.svg")

plt.show()
