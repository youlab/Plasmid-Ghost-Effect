import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
idx = [0, 4, 3]
colors = [colors[i] for i in idx] # indexing the colors to correspond to the strian 1, 15, and 76
def half_life(P0, b, k):
    a = b+k
    p0 = P0/100
    tau = np.log((2*a-b*p0)/(a-b*p0))/a
    return tau

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
    # 0. Remove rows with NaN in mu --------------------------------------
    nan_idx = np.isnan(mu)
    t   = t[~nan_idx]
    mu  = mu[~nan_idx]
    se  = se[~nan_idx]

    if len(t) < 2:
        print("not enough points left")
        return np.nan, np.nan                      # not enough points left

    # 1. Convert SD → SE and move to log space ---------------------------------
    g     = np.log(mu)
    se_g  = se / mu                               # delta: σ/μ
    g_star = np.log(mu_star)

    # 2. Locate the bracketing segment ----------------------------------------
    idx = np.where(g <= g_star)[0]                # first point below threshold
    if len(idx) == 0 or idx[0] == 0:
        print("threshold never reached")
        return np.nan, np.nan                     # threshold never reached

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

time = np.array([0, 3, 7, 12, 20, 28])

Strains     = [1, 15, 76]
inits = np.array([100, 90, 50])
fig,axes=plt.subplots(3, 3, figsize=(5, 5))
fig2, ax2 = plt.subplots(1, 1, figsize=(4, 1.9))
for i,strain in enumerate(Strains):
    abundance = np.load(f"./processed_data/{strain}+R388_mean.npy")
    se_abundance = np.load(f"./processed_data/{strain}+R388_std.npy")

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_ic = abundance.shape[0] # number of initial conditions
    bio_rep = abundance.shape[1] # number of biological replicates

    for j in range(n_ic):
        ax=axes[i,j]
        Mean=abundance[j,:,:]
        Se=se_abundance[j,:,:]
        for k in range(bio_rep):
            yk=Mean[k,:]
            ykerr=Se[k,:]
            non_nan = ~np.isnan(yk)
            x = time[non_nan]
            yk=yk[non_nan]
            ykerr=ykerr[non_nan]
            ax.plot(x,yk, color=colors[i], linewidth=1.5)
            ax.scatter(x, yk, s=20, color=colors[i])
            ax.errorbar(x,yk,ykerr,marker='None',linewidth=0,elinewidth=1,capsize=2,color="k")
        ax.hlines(abundance[j,0,0]/2, -1, 31, color='gray', linestyle='--', linewidth=1)
        ax.set_xlim([-1, 31])
        ax.set_xticks([0, 10, 20, 30])
        ax.set_ylim([5, 120])
        ax.set_yticks([10, 100])

        ax.set_yscale("log")
        if not (i==2 and j==0):
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xlabel("time (days)")
            ax.set_ylabel("P%",rotation=0,va="center",ha="right")
        if i==0:
            if j==0:
                ax.set_title("P$_0$=%i%%"%(inits[j]))
            else:
                ax.set_title("%i%%" % (inits[j]))
        if j==0:
            ax.text(x=0.1,y=0.1,s=f"#{strain}+R388",fontsize=13,transform=ax.transAxes)

    HL=np.zeros((n_ic,bio_rep))
    SE_HL=np.zeros((n_ic,bio_rep))
    persisting = False
    for j in range(n_ic):
        for k in range(bio_rep):
            t,t_se=loglinear_crossing_time(time,abundance[j,k,:],se_abundance[j,k,:],abundance[j,k,0]/2)
            HL[j,k]=t
            SE_HL[j,k]=t_se
    mask = np.isnan(HL)
    HL[mask] = np.max(time)
    median = np.median(HL,axis=1)

    print(f"{strain}+R388 half-lives (days):")
    for k in range(n_ic):
        mean_HL = np.mean(HL[k,:])
        sem_HL = np.std(HL[k,:], ddof=1) / np.sqrt(bio_rep)
        print(f"P0={inits[k]}%: {mean_HL:.2f} ± {sem_HL:.2f}")

    for k in range(n_ic):
        mean_HL = np.mean(HL[k,:])
        SD_HL = np.std(HL[k,:], ddof=1)
        print(f"{strain}+R388 half-life at P0% = {inits[k]}: {mean_HL:.1f}±{SD_HL:.1f}, n = 3")
        for j in range(bio_rep):
            ax2.scatter(inits[k],HL[k,j],s=100,color=colors[i],linewidth=1,edgecolors='black',zorder=k)

            ax2.errorbar(inits[k], HL[k,j], SE_HL[k,j], marker='o', markersize=0, elinewidth=1, linewidth=0,
                         capsize=2, color="k", zorder=k)
    ax2.plot(inits,median,c='k',lw=1.5,zorder=-10)

fig.tight_layout()
fig.savefig(f"./figures/KV2_Plasmids_Dynamics.png",dpi=300)
fig.savefig(f"./figures/KV2_Plasmids_Dynamics.svg",dpi=300)

for i, strain in enumerate(Strains):
    ax2.scatter([],[],s=100,color=colors[i],linewidth=1,edgecolors='black',label=f'#{strain}')
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),title='host strain',
           handlelength = 0.7)

ax2.set_xlim([40,105])
ax2.set_xticks([50, 90, 100])
ax2.set_ylim([0,31])
ax2.set_yticks([0, 15, 30])
ax2.set_yticklabels([0, 15, "> 30"])
ax2.set_ylabel(r"$\tau_{1/2}$ (days)")
ax2.set_xlabel("P$_0$%")
fig2.subplots_adjust(
    left=0.25,
    right=0.60875,
    bottom=0.23,
    top=0.95
)
fig2.savefig(f"./figures/KV2_halflife.png",dpi=300)
fig2.savefig(f"./figures/KV2_halflife.svg")
plt.show()