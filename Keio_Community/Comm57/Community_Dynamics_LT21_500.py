import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

dilution = 500
Plasmids = ["R388","RP4","pCU1","R6K"]
Abs = ["Trim","Tet","Carb","Strp"]
time=np.array([0,2,3,7,10,15,20])
Conditions=["NoAb","Ab",]
Labels=["no pulse","+Trim",]
Reps=3
n_show = 12
cmap = plt.get_cmap('Set3')
colors = [cmap(i) for i in np.linspace(0, 1, n_show)] + ["#B3B3B3"]

sum = 0
for plasmid in Plasmids:
    all_data = []
    for i in range(len(Conditions)):
        for j in range(Reps):
            df = pd.read_excel(f"./processed_data/organized_NGS/{plasmid}_{dilution}_dilution.xlsx",
                               sheet_name=f"{Conditions[i]}_Rep{j+1}",index_col=0)
            df=np.array(df)[:-2]#exclude rows "total" and "unmatched
            comp = df/np.sum(df,axis=0)
            all_data.append(comp)
    all_data = np.array(all_data)
    np.save(f"../composition_all/Comm57_{plasmid}_{dilution}.npy",all_data)
    sum += np.sum(all_data,axis=(0,2))
rank_idx = np.argsort(-sum)
donors = np.arange(4)
rank_idx = np.r_[donors, rank_idx[~np.isin(rank_idx, donors)]] # move the donor strains (1-4) to the front
strain_ids = pd.read_excel(f"./processed_data/organized_NGS/{plasmid}_{dilution}_dilution.xlsx",
                           sheet_name=f"NoAb_Rep1")["Barcodes"]
strain_ids = np.array(strain_ids)[:-2]
strain_ids = strain_ids[rank_idx]
print(strain_ids)
fig, axes = plt.subplots(2, 4, figsize=(10, 4))
for p,plasmid in enumerate(Plasmids):
    for i in range(len(Conditions)):
        ax=axes[i,p]
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
        ax.set_xticks(np.arange(0,21,10))
        ax.set_xlim([0,20])
        label = "no pulse" if i == 0 else f"+{Abs[p]}"
        ax.set_ylabel(label)
        if p==0 and i==0:
            ax.set_title(f"Comm57+{plasmid}")
        if p!=0 and i==0:
            ax.set_title(f"+{plasmid}")
        if p==0 and i==1:
            ax.set_xlabel("time (days)")
        else:
            ax.set_xticklabels([])
fig.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.9)
labels = list(strain_ids[0:n_show])+["Other"]
legend_handles = [
    Patch(facecolor=colors[i], label=labels[i])
    for i in range(n_show+1)
]

fig.legend(handles=legend_handles,title="strain id",loc="center right",bbox_to_anchor=(1,0.5))
fig.savefig(f"./figures/comm_dynamics_{dilution}dilution_all.png",dpi=300)

plt.show()
