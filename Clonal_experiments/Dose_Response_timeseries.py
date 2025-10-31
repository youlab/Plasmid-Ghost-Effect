import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

host = "MG1655"
plasmids=["pSC101","colE1","pUC"]
Ab=["Kan","Kan","Spect"]
markers=["o","^","s"]
ab_conditions=[0,2,7] # conditions for optional visualization
time = np.arange(0,21,1)
# ------------------------------------------------------------------
# 1. time series ------------------------------------------
# ------------------------------------------------------------------

for i,plasmid in enumerate(plasmids):
    fig1, ax1 = plt.subplots(1, 1, figsize=(2.05, 1.9))
    abundance = np.load(f"./LT_data_py/{plasmid}_dose_response_mean.npy")
    std_abundance = np.load(f"./LT_Data_py/{plasmid}_dose_response_std.npy")

    # define a hard floor for plasmid abundance (0.1%) for visualization in log scale
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_ab = abundance.shape[0] # number of different antibiotic pulse concentrations
    bio_rep = abundance.shape[1] # number of biological replicates

    mean_dyn = np.mean(abundance, axis=1)
    std_dyn = np.std(abundance, axis=1, ddof=1)

    for j,idx in enumerate(ab_conditions):
        Mean=mean_dyn[idx,:]
        Std=std_dyn[idx,:]

        ax1.plot(time,Mean, color=colors[i], linewidth=1.5, zorder=-10)
        ax1.scatter(time, Mean, s=80, color=colors[i], linewidth=1, edgecolors='black', marker=markers[j])
        ax1.errorbar(time,Mean,Std,marker='None',linewidth=0,elinewidth=1,capsize=2,color="k")
    ax1.fill_between(x=[2, 3], y1=[1, 1], y2=[150, 150], color="#808080", zorder=-20, alpha=0.3)
    ax1.set_xlim([-0.5, 21.5])
    ax1.set_xticks([0, 10,20])
    ax1.set_ylim([1, 150])
    ax1.set_yticks([1, 10, 100])
    ax1.set_yscale("log")

    ax1.set_title(plasmid,fontsize=15)
    if i==0:
        ax1.set_xlabel("time (days)")
        ax1.set_ylabel("P%", rotation=0, va="center", ha="right")

    fig1.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
    fig1.savefig(f"./figures/{plasmid}_dose_response_timeseries.png",dpi=300)
    fig1.savefig(f"./figures/{plasmid}_dose_response_timeseries.svg")

fig2,ax2=plt.subplots(1,1,figsize=(6,2))
labels = [0,"2$^{-5}$",r"1 ($\times$50 $\mu$g/mL)"]
for i in range(3):
    ax2.scatter([],[], s=80, color="#B3B3B3", linewidth=1, edgecolors='black', marker=markers[i], label=labels[i])
ax2.legend(ncol=3,handlelength=0.5)
fig2.savefig("./figures/dose_response_legend.svg")
plt.show()