import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

host = "MG1655"
plasmids=["pSC101","colE1","pUC"]
Ab=["Kan","Kan","Spect"]
time = np.arange(0,21,1)
# ------------------------------------------------------------------
# 1. time series ------------------------------------------
# ------------------------------------------------------------------
fig1, axes1 = plt.subplots(3, 8, figsize=(12, 5))
for i,plasmid in enumerate(plasmids):
    abundance = np.load(f"./LT_data_py/{plasmid}_dose_response_mean.npy")
    std_abundance = np.load(f"./LT_Data_py/{plasmid}_dose_response_std.npy")

    # define a hard floor for plasmid abundance (0.1%) for visualization in log scale
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_ab = abundance.shape[0] # number of different antibiotic pulse concentrations
    bio_rep = abundance.shape[1] # number of biological replicates

    for j in range(n_ab):
        ax1=axes1[i,j]
        Mean=abundance[j,:,:]
        Std=std_abundance[j,:,:]
        for k in range(bio_rep):
            yk=Mean[k,:]
            ykerr=Std[k,:]
            non_nan = ~np.isnan(yk)
            x = time[non_nan]
            yk=yk[non_nan]
            ykerr=ykerr[non_nan]
            ax1.plot(x,yk, color=colors[j], linewidth=1.5)
            ax1.errorbar(x,yk,ykerr,marker='None',linewidth=0,elinewidth=1,capsize=2,zorder=20,color="k")
            ax1.scatter(x,yk, s=20, color=colors[j], zorder=10)
        if j!=0:
            ax1.fill_between(x=[2, 3], y1=[1, 1], y2=[150, 150], color="#808080", zorder=-20, alpha=0.3)
        ax1.set_xlim([-1, 21])
        ax1.set_xticks([0, 10,20])
        ax1.set_ylim([1, 120])
        ax1.set_yticks([1, 10, 100])
        ax1.set_yscale("log")
        if not (i==2 and j==0):
            ax1.set_xticklabels([])
            ax1.set_yticklabels([])
        else:
            ax1.set_xlabel("time (days)")
            ax1.set_ylabel("P%",rotation=0,va="center",ha="right")
        bbox = dict(
            boxstyle="round,pad=0.3",  # round corners + padding
            facecolor="white",  # or any color
            edgecolor="none",  # or "gray", etc.
            alpha=0.8  # 0 = transparent, 1 = opaque
        )
        if j==0:
            s = plasmid+"\n"+r"[%s] = 0$\times$"%Ab[i]
            ax1.text(x=0.92,y=0.72,s=s,transform=ax1.transAxes,fontsize=11,ha="right",bbox=bbox,zorder=50)
        else:
            s = r"$2^{%s}\times$"%(-7+j)
            if j<3:
                ax1.text(x=0.92,y=0.72,s=s,transform=ax1.transAxes,fontsize=11,ha="right",bbox=bbox,zorder=50)
            else:
                ax1.text(x=0.12,y=0.12,s=s,transform=ax1.transAxes,fontsize=11,ha="left",bbox=bbox,zorder=50)

fig1.tight_layout()
fig1.savefig("./figures/dose_response_timeseries.png",dpi=300)
plt.show()