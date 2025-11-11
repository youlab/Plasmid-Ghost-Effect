import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.stats import t as t_stat

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
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

time = np.array([0, 3, 6, 7, 9, 15, 20, 25])

host = "DA28102"
plasmids=["pCU1","R6K"]
inits = np.array([100, 99, 90, 50, 20])
fig,axes=plt.subplots(2,5,figsize=(8,3.2))
for i,plasmid in enumerate(plasmids):
    abundance = np.load(f"./LT_data_py/{plasmid}_mean.npy")
    se_abundance = np.load(f"./LT_Data_py/{plasmid}_se.npy")

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
            ax.plot(x,yk, color=colors[j], linewidth=1.5)
            ax.scatter(x, yk, s=20, color=colors[j])
            ax.errorbar(x,yk,ykerr,marker='None',linewidth=0,elinewidth=1,capsize=2,color="k")
        ax.set_xlim([-1, 26])
        ax.set_xticks([0, 10, 20])
        ax.set_ylim([1, 120])
        ax.set_yticks([1, 10, 100])
        ax.set_yscale("log")
        if not (i==1 and j==0):
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
            ax.text(x=0.1,y=0.1,s=plasmid,fontsize=13,transform=ax.transAxes)
    fig2, ax2 = plt.subplots(1, 1, figsize=(2.05, 1.9))
    HL=np.zeros((n_ic,bio_rep))
    SE_HL=np.zeros((n_ic,bio_rep))
    persisting = False
    for j in range(n_ic):
        for k in range(bio_rep):
            t,t_se=loglinear_crossing_time(time,abundance[j,k,:],se_abundance[j,k,:],abundance[j,k,0]/2)
            HL[j,k]=t
            SE_HL[j,k]=t_se
    mask = np.isnan(HL)
    HL[mask] = 25
    median = np.median(HL,axis=1)

    for k in range(n_ic):
        mean_HL = np.mean(HL[k,:])
        SD_HL = np.std(HL[k,:], ddof=1)
        print(f"{plasmid} half-life at P0% = {inits[k]}: {mean_HL:.1f}±{SD_HL:.1f}, n = 3")
        for j in range(bio_rep):
            ax2.scatter(inits[k],HL[k,j],s=100,color=colors[k],linewidth=1,edgecolors='black',zorder=k)

            ax2.errorbar(inits[k], HL[k,j], SE_HL[k,j], marker='o', markersize=0, elinewidth=1, linewidth=0,
                         capsize=2, color="k",zorder=k)
    ax2.plot(inits,median,c='k',lw=1.5,zorder=-10)

    ax2.text(x=0.1, y=0.9, s=f"DA28102\n{plasmid}", fontsize=13, ha="left", va="top", transform=ax2.transAxes)
    ax2.set_xlim([10,105])
    ax2.set_xticks([50,100])

    if plasmid == "pCU1":
        ax2.set_ylim([0, 13])
        ax2.set_yticks([0,6,12])

    elif plasmid == "R6K":
        ax2.set_yticks([0,12,25])
        ax2.set_yticklabels([0,12,"> 25"])
        ax2.set_ylim([0,27])
    fig2.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
    fig2.savefig(f"./figures/{plasmid}_halflife.png",dpi=300)
    fig2.savefig(f"./figures/{plasmid}_halflife.svg")
fig.subplots_adjust(bottom=0.18)
fig.savefig("./figures/Conjugative_Plasmids_Plating.png",dpi=300)
plt.show()