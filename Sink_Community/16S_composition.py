import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Plasmids = ["pSC101","colE1","pUC"]
Abs = ["Kan","Kan","Spect"]
x = np.arange(4)

Conditions=["NoAb","Ab",]
Reps=3
cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
colors = [colors[0]] + colors[3:]

rank_idx = np.array([3,4,0,1,2,5])
genus_ids = pd.read_excel("./processed_data/genus_level_composition_filtered.xlsx",
                   sheet_name="Sheet1")["Genus"]
genus_ids = genus_ids.iloc[rank_idx].to_list()

for j,plasmid in enumerate(Plasmids):
    fig, ax = plt.subplots(1, 1, figsize=(2.05, 1.9))
    composition = np.load(f"./processed_data/{plasmid}_NoAb.npy")[:,rank_idx,:]
    Ab_comp = np.load(f"./processed_data/{plasmid}_Ab.npy")[:,rank_idx,-1]
    Ab_comp = Ab_comp[:,:,np.newaxis]
    composition = np.concatenate([composition,Ab_comp],axis=2)
    mean_composition = np.mean(composition,axis=0)
    std_composition = np.std(composition,axis=0,ddof=1)
    std_composition[:,0]=np.nan
    base = np.zeros(4)
    for k in range(5):
        mean = mean_composition[k,:]
        std = std_composition[k,:]
        se = std/np.sqrt(3)

        ax.errorbar(x, base+mean, se, marker='None', linewidth=0,
                    elinewidth=1, capsize=2, capthick=1, color="k")
        ax.bar(x, mean, bottom = base, color=colors[k])
        base = base+mean
    ax.set_yticks([])
    ax.set_ylim([0,1])
    ax.set_xticks(x)
    ax.set_title(plasmid)

    ax.set_xticklabels(["D0","D1","LB","+Ab"])

    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.87)
    fig.savefig(f"./figures/16S_{plasmid}.png",dpi=300)
    fig.savefig(f"./figures/16S_{plasmid}.svg")

from matplotlib.patches import Patch
fig2,ax2=plt.subplots(1,1)
labels = genus_ids+["Other"]
legend_handles = [
    Patch(facecolor=colors[i], label=labels[i])
    for i in range(5)
]

ax2.legend(handles=legend_handles,title="genus")
fig2.savefig(f"./figures/16S_legends.svg")
plt.show()
