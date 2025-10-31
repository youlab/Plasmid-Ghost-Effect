import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def p_to_star(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

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
    # 0. Remove rows with NaN in mu or sd --------------------------------------
    nan_idx = np.isnan(mu)
    t   = t[~nan_idx]
    mu  = mu[~nan_idx]
    se  = se[~nan_idx]

    if len(t) < 2:
        #print("not enough points left")
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

time = np.array([0, 1, 2, 3, 4, 5, 10])

Chemicals     = ["Rif", "Rif+Prm", "Rif+Mpt", "Rif+Phe"]
Chem_free_mean = np.load(f"./processed_data/NoChem_mean.npy")
Chem_free_se = np.load(f"./processed_data/NoChem_se.npy")
Treatment_length = ["1-day","2-day","3-day","continuous"]
fig,axes=plt.subplots(4,4,figsize=(6,6))
fig2, axes2 = plt.subplots(4, 1, figsize=(4, 6))
rng = np.random.default_rng(42)  # set seed for reproducible jitter
for i,chem in enumerate(Chemicals):
    ax2=axes2[i]
    abundance = np.load(f"./processed_data/{chem}_mean.npy")
    se_abundance = np.load(f"./processed_data/{chem}_se.npy")

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_tc = abundance.shape[0] # number of treatment duration conditions
    bio_rep = abundance.shape[1] # number of biological replicates

    for j in range(n_tc):
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
        # plot the chemical-free control
        for k in range(bio_rep):
            yk = Chem_free_mean[k, :]
            ykerr = Chem_free_se[k, :]
            non_nan = ~np.isnan(yk)
            x = time[non_nan]
            yk = yk[non_nan]
            ykerr = ykerr[non_nan]
            ax.plot(x, yk, color='#808080', linewidth=1.5)
            ax.scatter(x, yk, s=20, color='#808080')
            ax.errorbar(x, yk, ykerr, marker='None', linewidth=0, elinewidth=1, capsize=2, color='k')
        if j < 3:
            ax.fill_between(x=[0, j+1], y1=[1, 1], y2=[150, 150], color="#808080", zorder=-20, alpha=0.1)
        else:
            ax.fill_between(x=[0,10], y1=[1, 1], y2=[150, 150], color="#808080", zorder=-20, alpha=0.1)
        ax.set_xlim([-1, 11])
        ax.set_xticks([0, 5, 10])
        ax.set_ylim([10, 120])
        ax.set_yticks([10, 100])
        ax.set_yscale("log")
        if not (i==3 and j==0):
            ax.set_xticklabels([])
            ax.set_yticklabels([])
        else:
            ax.set_xlabel("time (days)")
            ax.set_ylabel("P%",rotation=0,va="center",ha="right")
        if i==0:
            ax.set_title(Treatment_length[j])
        if j==0:
            ax.text(x=0.1,y=0.1,s=chem,fontsize=13,transform=ax.transAxes)


    HL = np.zeros((n_tc+1, bio_rep))
    SE_HL = np.zeros((n_tc+1, bio_rep))
    for j in range(n_tc):
        Mean=abundance[j,:,:]
        Se=se_abundance[j,:,:]
        for k in range(bio_rep):
            t, t_se = loglinear_crossing_time(time, Mean[k, :], Se[k, :], 50)
            HL[j+1, k] = t
            SE_HL[j+1, k] = t_se
        for k in range(bio_rep):
            t, t_se = loglinear_crossing_time(time, Chem_free_mean[k, :], Chem_free_se[k, :], 50)
            HL[0, k] = t
            SE_HL[0, k] = t_se
    pvals = []
    pairs = [(0,4), (0,3), (0,2), (0,1)]
    for a,b in pairs:
        t_stat, p_value = stats.ttest_ind(HL[a, :], HL[b, :], equal_var=False)
        pvals.append(p_value)
    #
    # --- summary stats ---
    means = HL.mean(axis=1)
    sems = HL.std(axis=1, ddof=1) / np.sqrt(bio_rep)  # SEM across biological replicates
    print(f"Chemical-free group half-life: {means[0]:.1f}±{sems[0]:.1f} days; ")
    for j in range(4):
        print(f"{chem} treatment {Treatment_length[j]} half-life: {means[j+1]:.1f}±{sems[j+1]:.1f} days; ")
    # --- plot ---
    x = np.arange(5)

    # bars + error bars
    error_kw = {'elinewidth': 1.2, 'capsize': 3, 'capthick': 1.2, 'color': 'blue'}
    bars = ax2.bar(x, means, yerr=sems, width=0.6, error_kw=error_kw, lw=1.2, facecolor="None", edgecolor="k", zorder=2)

    # jittered points (strip-like)

    for j in range(5):
        jitter_scale = 0.1
        x_j = rng.normal(loc=x[j], scale=jitter_scale, size=bio_rep)
        color = "#808080" if j==0 else colors[i]
        ax2.scatter(x_j, HL[j], s=60, linewidths=1, facecolor=color, edgecolors='k', zorder=3)
        ax2.errorbar(x_j, HL[j], SE_HL[j], marker='o', markersize=0, elinewidth=1.2, linewidth=0,
                     capsize=2, capthick=1.2, color="k", zorder=4)
    for (a, b), p in zip(pairs, pvals):
        if p>=0.05:
            continue
        stars = p_to_star(p)
        y = 9.8
        ax2.text(b, y, stars, ha='center', va='bottom', fontsize=12)
    ax2.set_ylim([0, 13])
    ax2.set_yticks([0,5,10])
    ax2.set_xticks(x)
    ax2.set_title(chem)
    if i==3:
        ax2.set_xticklabels(["Chem-"]+Treatment_length,rotation=90)
        ax2.set_ylabel(r"$\tau_{1/2}$ (days)")
    else:
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])

fig.subplots_adjust(left=0.15, right=0.95, bottom=0.16, top=0.92, hspace=0.3)
fig.savefig("./figures/Chemical_Treatment_Plating.png",dpi=300)
fig.savefig("./figures/Chemical_Treatment_Plating.svg")
fig2.subplots_adjust(left=0.15, right=0.95, bottom=0.16, top=0.92, hspace=0.3)
fig2.savefig(f"./figures/treatment_half_lives.png", dpi=300)
fig2.savefig(f"./figures/treatment_half_lives.svg")
plt.show()