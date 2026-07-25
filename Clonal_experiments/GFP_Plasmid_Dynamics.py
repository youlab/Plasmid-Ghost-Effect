import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    """
    # 0. Remove rows with NaN in mu or se --------------------------------------
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

time = np.arange(0,13,1)

host = "MG1655"
plasmids=["pSC101","colE1","pUC"]
inits = np.array([100, 99, 90, 50, 20])
PCN_file = np.genfromtxt(
    "../qPCR/processed_data/PCN_summary_report.tsv",
    delimiter="\t",
    names=True,
    dtype=None,
    encoding="utf-8",
)
PCN_rows = np.atleast_1d(PCN_file)
PCN = {row["plasmid"]: float(row["mean"]) for row in PCN_rows}
SE_PCN = {row["plasmid"]: float(row["sem"]) for row in PCN_rows}

fig,axes=plt.subplots(3,5,figsize=(8,5))
for i,plasmid in enumerate(plasmids):
    abundance = np.load(f"./LT_data_py/{plasmid}_mean.npy")
    std_abundance = np.load(f"./LT_Data_py/{plasmid}_std.npy")

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_ic = abundance.shape[0] # number of initial conditions
    bio_rep = abundance.shape[1] # number of biological replicates

    for j in range(n_ic):
        ax=axes[i,j]
        Mean=abundance[j,:,:]
        Std=std_abundance[j,:,:]
        for k in range(bio_rep):
            yk=Mean[k,:]
            ykerr=Std[k,:]
            non_nan = ~np.isnan(yk)
            x = time[non_nan]
            yk=yk[non_nan]
            ykerr=ykerr[non_nan]
            ax.plot(x,yk, color=colors[j], linewidth=1.5)
            ax.scatter(x, yk, s=20, color=colors[j])
            ax.errorbar(x,yk,ykerr,marker='None',linewidth=0,elinewidth=1,capsize=2,color="k")

        ax.set_xlim([-1, 13])
        ax.set_xticks([0, 6, 12])
        ax.set_ylim([1, 120])
        ax.set_yticks([1, 10, 100])
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
            ax.text(x=0.1,y=0.1,s=plasmid,transform=ax.transAxes,fontsize=13)

    fig2, ax2 = plt.subplots(1, 1, figsize=(2.05, 1.9))
    HL=np.zeros((n_ic,bio_rep))
    SE_HL=np.zeros((n_ic,bio_rep))
    for j in range(n_ic):
        for k in range(bio_rep):
            t,t_se=loglinear_crossing_time(time,abundance[j,k,:],std_abundance[j,k,:],abundance[j,k,0]/2) # Here for a single replicate, std = se
            HL[j,k]=t
            SE_HL[j,k]=t_se

    mask = np.isnan(HL)
    HL[mask] = 12
    median = np.median(HL,axis=1)

    # save the per-replicate half-lives for each initial condition
    df = pd.DataFrame({
        "P0": np.repeat(inits[:n_ic], bio_rep),
        "replicate": np.tile(np.arange(1, bio_rep + 1), n_ic),
        "half_life": HL.ravel(),
        "half_life_se": SE_HL.ravel(),
        "right_censored": mask.ravel(),
    })
    df.to_csv(f"./LT_Data_py/{plasmid}_halflife.csv", index=False)

    print(f"{plasmid} half-lives (days):")
    for k in range(n_ic):
        print(f"{plasmid} half-life at P0% = {inits[k]}: ")
        print(HL[k,:])
        mean_HL = np.mean(HL[k,:])
        sd_HL = np.std(HL[k,:], ddof=1)
        sem_HL = sd_HL / np.sqrt(bio_rep)
        print(f"P0={inits[k]}%: {mean_HL:.2f} ± {sd_HL:.2f}, n = 3")
        if k==4:
            delta = np.log(2) / mean_HL
            se_delta = np.log(2) * sem_HL / (mean_HL**2)
            pcn = PCN[plasmid]
            se_pcn = SE_PCN[plasmid]
            kappa = 0.5**(pcn-1) * np.log2(500) / 24 # 500-fold daily dilution, converted to generation per hour
            max_HL = np.log(delta / kappa) / delta
            dmax_ddelta = (1.0 - np.log(delta / kappa)) / (delta**2)
            dmax_dpcn = np.log(2) / delta
            se_max_HL = np.sqrt((dmax_ddelta * se_delta)**2 + (dmax_dpcn * se_pcn)**2)
            print(
                f"fitness cost: {delta:.2f} ± {se_delta:.2f} day ^{-1}, "
                f"expected half-life at P0 = 100%: {max_HL:.3e} ± {se_max_HL:.3e} days"
            )


    for j in range(bio_rep):
        ax2.scatter(inits,HL[:,j],s=100,c=colors[0:5],zorder=10*bio_rep,linewidth=1,edgecolors='black')
        ax2.errorbar(inits, HL[:,j], SE_HL[:,j], marker='o', markersize=0, elinewidth=1, linewidth=0,
                     capsize=2, zorder=10*bio_rep, color="k")
    ax2.plot(inits,median,c='k',lw=1,zorder=-10)

    ax2.text(x=0.1,y=0.9,s=f"MG1655\n{plasmid}",fontsize=13,ha="left",va="top",transform=ax2.transAxes)

    ax2.set_xlim([10,105])
    ax2.set_xticks([50,100])

    if plasmid == "pSC101":
        ax2.set_ylim([0, 6.6])
        ax2.set_yticks([0, 3, 6])
        ax2.set_ylabel(r"$\tau_{1/2}$ (days)")
        ax2.set_xlabel("P$_0$%")
    elif plasmid == "colE1":
        ax2.set_ylim([0, 12])
        ax2.set_yticks([0,5,10])
    elif plasmid == "pUC":
        ax2.set_ylim([0, 13])
        ax2.set_yticks([0, 6, 12])
        ax2.set_yticklabels([0,6,"> 12"])

    fig2.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
    # fig2.savefig(f"./figures/{plasmid}_halflife.png",dpi=300)
    # fig2.savefig(f"./figures/{plasmid}_halflife.svg")
# fig.savefig("./figures/GFP_Plasmids_dynamics.png",dpi=300)
# fig.savefig("./figures/GFP_Plasmids_dynamics.svg")
plt.show()