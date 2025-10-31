import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

def loglinear_crossing_time(
        t: np.ndarray,
        mu: np.ndarray,
        se: np.ndarray
    ):
    """
    Estimate the time when mu(t) falls to `mu[0]/2` using log-linear
    interpolation and return its delta-method standard error (SE).

    Propagate the uncertainty in the plasmid abundance after the pulse (mu[0])

    Parameters
    ----------
    t   : 1-D array of time points (ascending).
    mu  : 1-D array of sample means at each time.
    se  : 1-D array of sample SDs (since there's 1 replicate, SD == SE).
    n_repl : int
        Number of raw replicates pooled in each (mu, sd).

    Returns
    -------
    t_cross : float
        Estimated crossing time.
    se_cross : float
        Standard error of `t_cross`.
        Returns (np.nan, np.nan) if the threshold is never crossed.
    """

    # ------------------------------------------------------------------ 1
    # SD → SE and log-transform -----------------------------------------
    g      = np.log(mu)
    se_g   = se / mu                        # delta method: σ/μ

    # Threshold (derived from first point) ------------------------------
    g0   = g[0]
    var_g0 = se_g[0] ** 2                   # Var[ln(mu[0])]
    g_star = np.log(0.5) + g0            # ln(factor * mu[0])
    var_gstar = var_g0                      # same variance

    # ------------------------------------------------------------------ 2
    # Locate bracketing segment -----------------------------------------
    idx = np.where(g <= g_star)[0]          # first point below threshold
    if idx.size == 0 or idx[0] == 0:
        #print("threshold never reached")
        return np.nan, np.nan

    i   = idx[0] - 1                        # segment [i, i+1]
    dt  = t[i + 1] - t[i]
    g1, g2 = g[i], g[i + 1]
    D   = g1 - g2
    N   = g1 - g_star

    t_cross = t[i] + dt * N / D             # log-linear interpolation

    # ------------------------------------------------------------------ 3
    # Delta-method variance ---------------------------------------------
    d_t_dg1   =  dt * (g_star - g2) / D**2
    d_t_dg2   =  dt * (g1     - g_star) / D**2
    d_t_dgstr = -dt / D                     # derivative wrt g_star

    var_t = (
        d_t_dg1**2   * se_g[i]**2 +
        d_t_dg2**2   * se_g[i + 1]**2 +
        d_t_dgstr**2 * var_gstar
    )

    se_cross = np.sqrt(var_t)
    return t_cross, se_cross


host = "MG1655"
plasmids=["pSC101","colE1","pUC"]
Ab=["Kan","Kan","Spect"]
time = np.arange(3,21,1)

ylims = [[0, 7.5], [0, 19], [0, 19]]
yticks = [[0, 3, 6], [0, 8, 17], [0, 8, 17]]
yticklabels = [[0, 3, 6], [0, 8, ">17"], [0, 8, ">17"]]
fold_dilution = np.array([1 / 128, 1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1]) # first element is zero, but here use 1/128 for visualization purpose
dilution_label = [0] + ["2$^{%i}$" % ii for ii in range(-6, 1, 2)]
# ------------------------------------------------------------------
# 2. half-lives ------------------------------------------
# ------------------------------------------------------------------
for i,plasmid in enumerate(plasmids):
    abundance = np.load(f"./LT_data_py/{plasmid}_dose_response_mean.npy")[:,:,3:]#ignore the time before the pulse
    std_abundance = np.load(f"./LT_Data_py/{plasmid}_dose_response_std.npy")[:,:,3:]#ignore the time before the pulse

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_ab = abundance.shape[0] # number of different antibiotic pulse concentrations
    bio_rep = abundance.shape[1] # number of biological replicates


    HL=np.zeros((n_ab,bio_rep))
    SE_HL=np.zeros((n_ab,bio_rep))
    for j in range(n_ab):
        for k in range(bio_rep):
            t,t_se=loglinear_crossing_time(time,abundance[j,k,:],std_abundance[j,k,:])
            HL[j,k]=t-3 # start counting after the antibiotic pulse
            SE_HL[j,k]=t_se

    fig1, ax1 = plt.subplots(1, 1, figsize=(2.05, 1.9))
    HL[np.isnan(HL)]=17
    for j in range(bio_rep):
        tau = HL[:,j]
        se = SE_HL[:,j]
        ax1.scatter(fold_dilution, tau, s=80, color=colors[i], linewidth=1, edgecolors='black')
        ax1.errorbar(fold_dilution, tau, se, marker='o', markersize=0, elinewidth=1.5,linewidth=0,
                     capsize=2, color="k")

    median = np.median(HL,axis=1)
    ax1.plot(fold_dilution,median,lw=1.5,c=colors[i],zorder=-10)
    #ax1.text(x=0.98,y=0.05,s=plasmid,ha="right",va="bottom",transform=ax1.transAxes,fontsize=15)
    ax1.set_xscale("log", base=2)
    ax1.set_xticks([1 / 128, 1 / 64, 1 / 16, 1 / 4, 1])
    ax1.set_ylabel(r"$\tau_{1/2}$ (days)")
    ax1.set_xticklabels(dilution_label)
    ax1.set_xlabel(Ab[i] + r" ($\times$50 $\mu$g/mL)")
    ax1.set_ylim(ylims[i])
    ax1.set_yticks(yticks[i])
    ax1.set_yticklabels(yticklabels[i])
    ax1.spines[['right', 'top']].set_visible(False)
    fig1.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
    fig1.savefig(f"./figures/dose_response_{plasmid}.png", dpi=300)
    fig1.savefig(f"./figures/dose_response_{plasmid}.svg")
plt.show()