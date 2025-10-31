import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t

# Linear model
def linear_model(x, m, b):
    return m * x + b

def mean_and_var_of_mean(arr):
    """Return (mean, variance_of_the_mean) for a 1-D vector."""
    n   = arr.size
    var = np.var(arr, ddof=1) / n      # divide by n → variance *of the mean*
    return arr.mean(), var

# calibration based on GFP and OD reads at different mix ratio of plasmid-carrying and plasmid-free cells
mix_ratio = np.array([100,99,90,50,20,0])

# Focusing only on the calibration curve without antibiotic treatment
# ------------------------------------------------------------
# 1.  Load raw plate data
# ------------------------------------------------------------
od_df  = pd.read_excel('./Raw_data/GE_LT5_GFP.xlsx',
                       sheet_name='Calibration',
                       header=None, skiprows=2,  nrows=6, usecols="B:G")
gfp_df = pd.read_excel('./Raw_data/GE_LT5_GFP.xlsx',
                       sheet_name='Calibration',
                       header=None, skiprows=17, nrows=6, usecols="B:G")

od_np  = od_df.to_numpy()
gfp_np = gfp_df.to_numpy()

# ------------------------------------------------------------
# 2.  Blank statistics (3 blank wells)
# ------------------------------------------------------------
od_blank   = od_np [-1, 3:6]
gfp_blank  = gfp_np[-1, 3:6]

mean_od_blank,  var_od_blank  = mean_and_var_of_mean(od_blank)
mean_gfp_blank, var_gfp_blank = mean_and_var_of_mean(gfp_blank)

# instrument/well noise (= SD of a single well, used later)
var_od_well  = np.var(od_blank,  ddof=1)          # NOT divided by 3
var_gfp_well = np.var(gfp_blank, ddof=1)

# ------------------------------------------------------------
# 3.  MG1655 (plasmid-free control)
# ------------------------------------------------------------
MG_od        = od_np [-1, 0:3]
MG_gfp       = gfp_np[-1, 0:3]

mean_MG_od,  var_MG_od  = mean_and_var_of_mean(MG_od)
mean_MG_gfp, var_MG_gfp = mean_and_var_of_mean(MG_gfp)

# subtract blank means  (variance of a difference of two means = sum)
mean_MG_od  -= mean_od_blank
var_MG_od   += var_od_blank
std_MG_od    = np.sqrt(var_MG_od)

mean_MG_gfp -= mean_gfp_blank
var_MG_gfp  += var_gfp_blank
std_MG_gfp   = np.sqrt(var_MG_gfp)

# GFP ⁄ OD for MG1655
y_MG     = mean_MG_gfp / mean_MG_od
std_MG_y = y_MG * np.sqrt((std_MG_gfp / mean_MG_gfp) ** 2 +
                          (std_MG_od  / mean_MG_od ) ** 2)

# ------------------------------------------------------------
# 4.  Calibration wells (P0 = 20 - 100)
# ------------------------------------------------------------
od_corr  = od_np [0:5, 0:3] - mean_od_blank
gfp_corr = gfp_np[0:5, 0:3] - mean_gfp_blank

# variance of each well’s corrected reading = instrument + blank-mean part
var_od_tot  = var_od_well  + var_od_blank
var_gfp_tot = var_gfp_well + var_gfp_blank

# ratio GFP/OD and its SD for every well
y       = gfp_corr / od_corr
std_y   = np.sqrt(var_gfp_tot / od_corr**2 +
                  (gfp_corr**2) * var_od_tot / od_corr**4)

# Now, perform linear regression on each fo the three plasmids: pSC101, colE1, pUC
plasmids = ["pSC101","colE1","pUC"]
fit_result = {plasmid:[] for plasmid in plasmids}
fig,axes=plt.subplots(1,3,figsize=(8,3))
for i, plasmid in enumerate(plasmids):
    ax=axes[i]
    yp = np.append(y[:,i],y_MG)
    std_yp = np.append(std_y[:,i],std_MG_y)
    ax.errorbar(mix_ratio,yp,std_yp,marker="o",c='steelblue',mfc="None",capsize=2,lw=0,markersize=10,
                    elinewidth=1)
    # Weighted fit using standard deviation (since the input are single measurements)
    popt, pcov = curve_fit(linear_model, mix_ratio, yp, sigma=std_yp, absolute_sigma=True)
    slope, intercept = popt
    slope_se, intercept_se = np.sqrt(np.diag(pcov))

    fit_result[plasmid]=[slope,intercept,pcov]

    # Compute t-statistic and p-value
    dof = len(mix_ratio) - 2
    t_stat = slope / slope_se
    p_value = 2 * (1 - t.cdf(np.abs(t_stat), df=dof))
    # residuals and total sums of squares
    y_pred = linear_model(mix_ratio, slope, intercept)
    SS_res = np.sum(((yp - y_pred) / std_yp) ** 2)  # weighted residual SS
    SS_tot = np.sum(((yp - np.average(yp, weights=1 / std_yp ** 2)) / std_yp) ** 2)

    R2 = 1 - SS_res / SS_tot

    # --- 95 % confidence band -----------------------------------------------
    t_crit = t.ppf(0.975, dof)  # two-tailed 95 %

    x_band = np.linspace(-5, 105, 300)  # dense grid for smooth band
    y_band = linear_model(x_band, slope, intercept)

    # Jacobian of y wrt [slope, intercept] at every x -> [[x, 1], ...]
    J = np.vstack((x_band, np.ones_like(x_band))).T  # shape (N, 2)

    # point-wise variance:  diag( J · pcov · Jᵀ )
    sigma2_y = np.einsum('ij,jk,ik->i', J, pcov, J)
    ci_y = t_crit * np.sqrt(sigma2_y)  # half-width of CI

    # plot fit line and its confidence envelope
    ax.plot(x_band, y_band, c='k', lw=1)
    ax.fill_between(x_band, y_band - ci_y, y_band + ci_y,
                    color='steelblue', alpha=0.50, linewidth=0)

    # Format regression equation
    txt = f"\ny = {slope:.0f} x {'+' if intercept >= 0 else '-'} {abs(intercept):.0f}\nR$^2$ = {R2:.3f}"
    ax.text(x=0.1, y=0.78, s=txt, transform=ax.transAxes, fontsize=12)

    ax.set_xlim([-5, 105])
    ax.set_title(plasmid)
    # ----- formatted report ---------------------------------------------------
    print(f"{plasmid} calibration without antibiotic treatment")
    report = (
        f"Slope            : {slope:.3f} ± {slope_se:.3f} (SE)\n"
        f"Intercept        : {intercept:.3f} ± {intercept_se:.3f} (SE)\n"
        f"t-statistic      : {t_stat:.2f}  (dof = {dof})\n"
        f"p-value (slope≠0): {p_value:.3g}\n"
        f"Weighted R²      : {R2:.3f}"
    )

    print(report)
    print()
ax=axes[0]
ax.set_ylabel("GFP/OD")
ax.set_xlabel("P%")
fig.tight_layout()
fig.savefig("./figures/LT5_GFP_calibration.png",dpi=300)
fig.savefig("./figures/LT5_GFP_calibration.svg")

def calculate_abundance(y,var_y,plasmid):
    # m is the slope, and b is the intercept
    m,b,pcov = fit_result[plasmid]
    x = (y-b)/m
    var_m = pcov[0, 0]#slope
    var_b = pcov[1, 1]#intercept
    cov_mb = pcov[0, 1]#covariance

    sigma_x = np.sqrt(
        (var_y + var_b) / m ** 2 +
        (y - b) ** 2 * var_m / m ** 4 +
        2 * (y - b) * cov_mb / m ** 3
    ) # standard deviation of plasmid abundance

    return x, sigma_x

# Initialization
days = np.arange(1, 13)
bio_rep = 3
inits_condits = 5
plasmids = ["pSC101", "colE1", "pUC"]
T = len(days)+1
L = len(plasmids)

Mean_P = np.full((L, inits_condits, bio_rep, T), np.nan)
Std_P = np.full((L, inits_condits, bio_rep, T), np.nan)

# fill Day 0
for i,P0 in enumerate(mix_ratio[:-1]):
    Mean_P[:,i,:,0]=P0
    Std_P[:,i,:,0]=0

# Load data from Excel
for t in days:
    # load OD and GFP raw data
    od_df = pd.read_excel('./Raw_data/GE_LT5_GFP.xlsx', sheet_name=f'Day{t}', header=None, skiprows=2, nrows=6,
                          usecols="B:J")
    od_np = od_df.to_numpy()
    gfp_df = pd.read_excel('./Raw_data/GE_LT5_GFP.xlsx', sheet_name=f'Day{t}', header=None, skiprows=11, nrows=6,
                          usecols="B:J")
    gfp_np = gfp_df.to_numpy()

    # Blank statistics
    od_blank = od_np[-1, 3:6]
    gfp_blank = gfp_np[-1, 3:6]

    mean_od_blank, var_od_blank = mean_and_var_of_mean(od_blank)
    mean_gfp_blank, var_gfp_blank = mean_and_var_of_mean(gfp_blank)
    # instrument/well noise (= SD of a single well, used later)
    var_od_well = np.var(od_blank, ddof=1)  # NOT divided by 3
    var_gfp_well = np.var(gfp_blank, ddof=1)

    # reading after subtracting the blank
    od_corr = od_np[0:5, :] - mean_od_blank
    gfp_corr = gfp_np[0:5, :] - mean_gfp_blank

    # variance of each well’s reading = instrument + blank-mean part
    var_od_tot = var_od_well + var_od_blank
    var_gfp_tot = var_gfp_well + var_gfp_blank

    # ratio GFP/OD and its variance for every well
    y = gfp_corr / od_corr
    var_y = var_gfp_tot / od_corr ** 2 + (gfp_corr ** 2) * var_od_tot / od_corr ** 4

    for l, plasmid in enumerate(plasmids):
        start_col = l * bio_rep
        end_col = (l + 1) * bio_rep
        p, std_p = calculate_abundance(y[:,start_col:end_col], var_y[:,start_col:end_col], plasmid)# estimated plasmid abundance
        Mean_P[l, :, :, t] = p
        Std_P[l, :, :, t] = std_p
# ---------------------------------------------------------------------
# save one *_mean.npy and *_std.npy per plasmid
# ---------------------------------------------------------------------
for l, plasmid in enumerate(plasmids):
    abundance = Mean_P[l, :, :, :]         # shape (inits_condits, bio_rep, T)
    std       = Std_P[l, :, :, :]

    np.save(f"./LT_data_py/{plasmid}_mean.npy",abundance)
    np.save(f"./LT_data_py/{plasmid}_std.npy", std)

plt.show()
