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

plasmids=["pSC101","colE1","pUC"]

for i,plasmid in enumerate(plasmids):
    HL=[]
    SE_HL=[]
    # process LT5 dynamics first, focus only on the 50% group
    time = np.arange(0,13,1)
    abundance = np.load(f"./LT_data_py/{plasmid}_mean.npy")
    std_abundance = np.load(f"./LT_Data_py/{plasmid}_std.npy")

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    for k in range(3):
        t,t_se=loglinear_crossing_time(time,abundance[-2,k,:],std_abundance[-2,k,:],abundance[-2,k,0]/2) # Here for a single replicate, std = se
        HL.append(t)
        SE_HL.append(t_se)

    time = np.arange(0,21,1)
    abundance = np.load(f"./LT_data_py/{plasmid}_dose_response_mean.npy")
    std_abundance = np.load(f"./LT_Data_py/{plasmid}_dose_response_std.npy")

    # define a hard floor for plasmid abundance (0.1%), does not affect the estimation of half-lives, but as a quick fix
    # for performing log-linear interpolation
    zero_idx = (abundance<0.1)
    abundance[zero_idx]=0.1

    n_ab = abundance.shape[0] # number of different antibiotic pulse concentrations
    bio_rep = abundance.shape[1] # number of biological replicates
    for k in range(3):
        t,t_se=loglinear_crossing_time(time,abundance[0,k,:],std_abundance[0,k,:],abundance[0,k,0]/2) # Here for a single replicate, std = se
        HL.append(t)
        SE_HL.append(t_se)

    mean_HL = np.mean(HL)
    std_HL = np.std(HL, ddof=1)
    print(f"ancestral {plasmid} half-life at P0% = 50%: {mean_HL:1.1f} ± {std_HL:1.1f} days")