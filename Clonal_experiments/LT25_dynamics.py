import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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

evolved_time = 3 # 3 or 20
plasmids=["pSC101","colE1","pUC"]

df = {"plasmid":[],"half_life":[],"half_life_se":[]}
for i,plasmid in enumerate(plasmids):
    fig,axes=plt.subplots(3,3,figsize=(5, 5))
    axes=axes.flat
    abundance = np.load(f"./LT_Data_py/LT_25_{plasmid}_day{evolved_time}_mean.npy")
    std_abundance = np.load(f"./LT_Data_py/LT_25_{plasmid}_day{evolved_time}_std.npy")

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_c = abundance.shape[0] # number of clones
    bio_rep = abundance.shape[1] # number of biological replicates

    for j in range(n_c):
        ax=axes[j]
        Mean=abundance[j,:,:]
        Std=std_abundance[j,:,:]
        time = np.arange(0,Mean.shape[1],1)
        for k in range(bio_rep):
            yk=Mean[k,:]
            ykerr=Std[k,:]
            non_nan = ~np.isnan(yk)
            x = time[non_nan]
            yk=yk[non_nan]
            ykerr=ykerr[non_nan]
            ax.plot(x,yk, color=colors[i], linewidth=1.5)
            ax.scatter(x, yk, s=20, color=colors[i])
            ax.errorbar(x,yk,ykerr,marker='None',linewidth=0,elinewidth=1,capsize=2,color="k")
        ax.plot([-1, time[-1]+1], [25, 25], 'k--', lw=1)
        ax.set_xlim([-1, time[-1]+1])
        ax.set_xticks([0, 2, 4])
        ax.set_ylim([1, 120])
        ax.set_yticks([1, 10, 100])
        ax.set_yscale("log")
        if not j == 6:
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xlabel("time (days)")
            ax.set_ylabel(f"{plasmid}%")
    fig.tight_layout()
    fig.savefig(f"./figures/LT25_{plasmid}_day{evolved_time}_comp_assay.png",dpi=300)
    fig.savefig(f"./figures/LT25_{plasmid}_day{evolved_time}_comp_assay.svg")

    HL=np.zeros((n_c,bio_rep))
    SE_HL=np.zeros((n_c,bio_rep))
    for j in range(n_c):
        for k in range(bio_rep):
            t,t_se=loglinear_crossing_time(time,abundance[j,k,:],std_abundance[j,k,:],abundance[j,k,0]/2) # Here for a single replicate, std = se
            HL[j,k]=t
            SE_HL[j,k]=t_se
    
    HL[np.isnan(HL)]=time[-1]
    mean_HL = np.nanmean(HL, axis=1)
    std_HL = np.nanstd(HL, axis=1, ddof=1)
    # Set standard deviation to NaN for right-censored data (where half-life is the maximum observed time)
    idx = np.where((mean_HL == 4) & (std_HL == 0))[0]
    std_HL[idx] = np.nan

    df["plasmid"].extend([plasmid]*n_c)
    df["half_life"].extend(mean_HL)
    df["half_life_se"].extend(std_HL)
fig2,ax2 = plt.subplots(1,1,figsize=(3,1.8))
df = pd.DataFrame(df)
df.to_csv(f"./LT_Data_py/LT25_day{evolved_time}_halflife.csv", index=False)

order = plasmids
palette = dict(zip(order, colors[:len(order)]))

sns.boxplot(
    data=df,
    x="plasmid",
    y="half_life",
    order=order,
    ax=ax2,
    width=0.55,
    showfliers=False,
    boxprops={"facecolor": "none", "edgecolor": "k", "linewidth": 1.2},
    whiskerprops={"color": "k", "linewidth": 1.2},
    capprops={"color": "k", "linewidth": 1.2},
    medianprops={"color": "k", "linewidth": 1.2},
    zorder=1,
)

summary = (
    df.groupby("plasmid")["half_life"]
    .agg(
        median="median",
        q1=lambda x: x.quantile(0.25),
        q3=lambda x: x.quantile(0.75),
    )
)

summary["IQR"] = summary["q3"] - summary["q1"]

# keep same order as boxplot
summary = summary.loc[order]

print(f"plasmids evolved until day {evolved_time} half-life summary:")
print(summary.to_string(float_format="{:.1f}".format))

jitter_scale = 0.08
rng = np.random.default_rng(42)
for i, plasmid in enumerate(plasmids):
    df_i = df[df["plasmid"] == plasmid]
    x_i = rng.normal(loc=i, scale=jitter_scale, size=len(df_i))
    ax2.scatter(x_i, df_i["half_life"], s=100, linewidths=1, facecolor=colors[i], edgecolors='k', zorder=3)
    ax2.errorbar(x_i, df_i["half_life"], yerr=df_i["half_life_se"], marker='o', markersize=0, elinewidth=1.2, linewidth=0,
                    capsize=2, capthick=1.2, color="k", zorder=4)
ax2.set_ylabel(r"$\tau_{1/2}$ (days)")
ax2.set_yticks([0, 2, 4])
if evolved_time == 3:
    ax2.set_yticklabels([0, 2, ">4"])
ax2.set_ylim([0,4.5])
ax2.set_xlabel("")
ax2.set_title(f"day {evolved_time}")
fig2.tight_layout()
fig2.savefig(f"./figures/LT25_day{evolved_time}_halflife.png",dpi=300)
fig2.savefig(f"./figures/LT25_day{evolved_time}_halflife.svg")

plt.show()