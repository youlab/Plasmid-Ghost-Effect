import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

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

def loglinear_crossing_time(
        t: np.ndarray,
        mu: np.ndarray,
        se: np.ndarray,
        mu_star: float
    ):
    """
    Log-linear interpolation on ln(mu) with delta-method SE,
    robust to NaNs in mu or sd (rows with NaNs are skipped).
    """
    # 0. Remove rows with NaN in mu or sd --------------------------------------
    nan_idx = np.isnan(mu)
    t   = t[~nan_idx]
    mu  = mu[~nan_idx]
    se  = se[~nan_idx]

    if len(t) < 2:
        print("not enough points left")
        return np.nan, np.nan                      # not enough points left

    # 1. Move to log space ---------------------------------
    g     = np.log(mu)
    se_g  = se / mu                               # delta: σ/μ
    g_star = np.log(mu_star)

    # 2. Locate the bracketing segment ----------------------------------------
    idx = np.where(g <= g_star)[0]                # first point below threshold
    if len(idx) == 0 or idx[0] == 0:
        #print("threshold never reached")
        return t[-1], np.nan                     # threshold never reached

    i  = idx[0] - 1                               # segment [i, i+1]
    dt = t[i + 1] - t[i]
    N  = g[i] - g_star
    D  = g[i] - g[i + 1]

    t_cross = t[i] + dt * N / D                   # log-linear interpolation

    # 3. Delta-method SE -------------------------------------------------------
    se_i, se_ip1 = se_g[i], se_g[i + 1]
    var_t  = (dt**2 / D**4) * ((g_star - g[i + 1])**2 * se_i**2 +
                               N**2 * se_ip1**2)
    se_cross = np.sqrt(var_t)
    return t_cross, se_cross

Plasmids = ["pSC101","colE1","pUC"]
Abs = ["Kan","Kan","Spect"]
conditions = ["+","-"] # with / without sponge
for i, plasmid in enumerate(Plasmids):
    baselines = np.load(f"./processed_data/Ecoli+{plasmid}.npy")
    for j,condition in enumerate(conditions):
        abundance = np.load(f"./processed_data/{plasmid}_mean.npy")[2*j:2*(j+1)]
        std_abundance = np.load(f"./processed_data/{plasmid}_std.npy")[2*j:2*(j+1)]
        baseline_j = baselines[j,:]
        fig,ax=plt.subplots(1,1,figsize=(2.05, 1.9))
        # define a hard floor for GFP/OD (1)
        zero_idx = (abundance < 1)
        abundance[zero_idx] = 0.1
        time = np.arange(0,11,1)

        n_ab = abundance.shape[0]  # with / without antibiotic pulse
        bio_rep = abundance.shape[1]  # number of biological replicates

        for p in range(n_ab):
            zorder = -p
            for q in range(bio_rep):
                ax.plot(time, abundance[p, q, :], color=colors[p], linewidth=1.5, zorder=zorder)
                ax.scatter(time, abundance[p, q, :], s=50, color=colors[p], linewidth=1, edgecolors='black', zorder=zorder)
                ax.errorbar(time, abundance[p, q, :], std_abundance[p ,q, :], marker='None', linewidth=0, elinewidth=1,
                            capsize=1.3, color="k", zorder=zorder)
        # save the plotted time series of each biological replicate (index 0 = LB, index 1 = LB+Ab)
        records = {"time": time}
        for p in range(n_ab):
            cond = "noAb" if p == 0 else "Ab"
            for q in range(bio_rep):
                records[f"mean_{cond}_rep{q+1}"] = abundance[p, q, :]
                records[f"se_{cond}_rep{q+1}"] = std_abundance[p, q, :]
        pd.DataFrame(records).to_csv(f"./processed_data/{plasmid}_{condition}sponge_dynamics.csv", index=False)

        ax.fill_between(x=[1, 2], y1=[1, 1], y2=np.ones(2)*baseline_j[0]**(1/0.9), color="#808080", zorder=-20, alpha=0.3)
        ax.plot([-1,11],np.ones(2)*baseline_j[0],c="#7FBF7B",lw=1,zorder=-9)
        ax.fill_between(x=[-1,11],y1=baseline_j[0]-baseline_j[1],y2=baseline_j[0]+baseline_j[1],color="#7FBF7B", zorder=-10, alpha=0.3)
        ax.text(x=0.97,y=0.97,s=plasmid,c="#7FBF7B",ha="right",va="top",transform=ax.transAxes)
        #ax.text(x=0.94, y=0.1, s=f"{plasmid}\n{condition}sponge", ha="right", transform=ax.transAxes)
        ax.set_xlim([-1, 11])
        ax.set_xticks([0, 5, 10])
        ax.set_yscale("log")
        ax.set_ylim([1e2,baseline_j[0]**(1/0.9)])
        if plasmid == "pSC101":
            ax.set_xlabel("time (days)")
            ax.set_ylabel("GFP/OD")
        fig.subplots_adjust(left=0.25, right=0.95, bottom=0.23, top=0.95)
        # fig.savefig(f"./figures/{plasmid}_{condition}sponge.png",dpi=300)
        # fig.savefig(f"./figures/{plasmid}_{condition}sponge.svg")
plt.show()
