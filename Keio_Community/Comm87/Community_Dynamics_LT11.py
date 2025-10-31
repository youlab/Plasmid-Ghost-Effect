import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plasmid = "R388"
Ab = "Trim"
dilution = 100

time=np.array([0,2,3,7,10,15,20,25,30,35])
Conditions=["NoAb","Ab",]
Labels=["no pulse","+Trim",]
Reps=3
cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
colors = [colors[0]] + colors[3:]

all_data = []
for i in range(len(Conditions)):
    for j in range(Reps):
        df = pd.read_excel("./processed_data/R388_100_dilution_exclude_zeros.xlsx",
                           sheet_name=f"{Conditions[i]}_Rep{j+1}",index_col=0)
        df=np.array(df)
        comp = df/np.sum(df,axis=0)
        all_data.append(comp)
all_data = np.array(all_data)
np.save("../composition_all/Comm87.npy",all_data)
sum = np.sum(all_data,axis=(0,2))
rank_idx = np.argsort(-sum)
rank_idx = np.r_[0, rank_idx[rank_idx != 0]] # move the donor strain to the front
strain_ids = pd.read_excel("./processed_data/R388_100_dilution_exclude_zeros.xlsx",
                   sheet_name="NoAb_Rep1")["Barcodes"]
strain_ids = np.array(strain_ids)[rank_idx]

for i in range(len(Conditions)):
    fig, ax = plt.subplots(1, 1, figsize=(2.05, 1.9))
    composition = []
    for j in range(Reps):
        df = pd.read_excel("./processed_data/R388_100_dilution_exclude_zeros.xlsx",
                           sheet_name=f"{Conditions[i]}_Rep{j+1}",index_col=0)
        df=np.array(df)
        comp = df/np.sum(df,axis=0)
        composition.append(comp)
    composition = np.array(composition)
    composition = composition[:,rank_idx,:]
    mean_composition = np.mean(composition,axis=0)
    std_composition = np.std(composition,axis=0,ddof=1)
    n = mean_composition.shape[0]
    base = np.zeros(len(time))
    for k in range(5):
        mean = mean_composition[k,:]
        std = std_composition[k,:]
        se = std/np.sqrt(3)
        ax.errorbar(time, base+mean, se, marker='None', linewidth=0.5,
                    elinewidth=0.5, capsize=2, capthick=0.5, color="k")
        ax.fill_between(time, base, base+mean, color=colors[k])
        base = base+mean
    top5_comp = np.sum(mean_composition[0:5,:],axis=0)
    ax.fill_between(time, top5_comp, np.ones_like(time), color=colors[-1])
    ax.set_yticks([])
    ax.set_ylim([0,1])
    ax.set_xticks(np.arange(0,35,10))
    ax.set_xlim([0,35])
    ax.set_title(Labels[i])

    if i==0:
        ax.set_ylabel("Comm87")
        ax.set_xlabel("time (days)")
    fig.subplots_adjust(left=0.25, right=0.95, bottom=0.15, top=0.87)
    fig.savefig(f"./figures/LT11_comm_dynamics_{Conditions[i]}.png",dpi=300)
    fig.savefig(f"./figures/LT11_comm_dynamics_{Conditions[i]}.svg")

from matplotlib.patches import Patch
fig2,ax2=plt.subplots(1,1)
labels = list(strain_ids[0:5])+["Other"]
legend_handles = [
    Patch(facecolor=colors[i], label=labels[i])
    for i in range(6)
]

ax2.legend(handles=legend_handles,title="strain id")
fig2.savefig(f"./figures/LT11_legends.svg")

plt.show()
