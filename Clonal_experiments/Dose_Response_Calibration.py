import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from scipy.stats import t
from sklearn.metrics import r2_score
# Linear model
def linear_model(x, m, b):
    return m * x + b

def mean_and_var_of_mean(arr):
    """Return (mean, variance_of_the_mean) for a 1-D vector."""
    n   = arr.size
    var = np.var(arr, ddof=1) / n      # divide by n → variance *of the mean*
    return arr.mean(), var

# calibration based on GFP and OD reads at different mix ratio of plasmid-carrying and plasmid-free cells
mix_ratio = np.array([100,80,50,20,0])

# ------------------------------------------------------------
# 1.  Load raw plate data
# ------------------------------------------------------------
od_df  = pd.read_excel('./Raw_data/GE_LT18_GFP_calibration.xlsx',
                       sheet_name='Calibration',
                       header=None, skiprows=2,  nrows=6, usecols="B:G")
gfp_df = pd.read_excel('./Raw_data/GE_LT18_GFP_calibration.xlsx',
                       sheet_name='Calibration',
                       header=None, skiprows=17, nrows=6, usecols="B:G")

od_np  = od_df.to_numpy()
gfp_np = gfp_df.to_numpy()

selection = "-" # +: with selection; -: without selection
start_col = 0 if selection =="+" else 3

# Calibration curve with antibiotic treatment
# ------------------------------------------------------------
# 2.  Blank and MG1655 (plasmid-free control)
# ------------------------------------------------------------
od_blank   = od_np[-1, 0]
gfp_blank  = gfp_np[-1, 0]

MG_od      = od_np [4, start_col] - od_blank
MG_gfp     = gfp_np[4, start_col] - gfp_blank
y_MG = MG_gfp/MG_od

# ------------------------------------------------------------
# 3.  Calibration wells (P0 = 20 - 100)
# ------------------------------------------------------------
od_corr  = od_np [0:4, start_col:start_col+3] - od_blank
gfp_corr = gfp_np[0:4, start_col:start_col+3] - gfp_blank

# ratio GFP/OD
y       = gfp_corr / od_corr

# Now, perform linear regression on each fo the three plasmids: pSC101, colE1, pUC
plasmids = ["pSC101","colE1","pUC"]
fit_result = {plasmid:[] for plasmid in plasmids}
fig,axes=plt.subplots(1,3,figsize=(8,3.3))
for i, plasmid in enumerate(plasmids):
    ax=axes[i]
    yp = np.append(y[:,i],y_MG)
    ax.scatter(mix_ratio,yp,marker="o",ec='steelblue',fc="None",s=50,lw=1.5)
    # Weighted fit using standard deviation (since the input are single measurements)
    popt, pcov = curve_fit(linear_model, mix_ratio, yp)
    slope, intercept = popt
    slope_se, intercept_se = np.sqrt(np.diag(pcov))

    fit_result[plasmid]=[slope,intercept,pcov]

    # Compute t-statistic and p-value
    dof = len(mix_ratio) - 2
    t_stat = slope / slope_se
    p_value = 2 * (1 - t.cdf(np.abs(t_stat), df=dof))
    # residuals and total sums of squares
    y_pred = linear_model(mix_ratio, slope, intercept)

    R2 = r2_score(yp,y_pred)

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
    print(f"{plasmid} calibration {selection} antibiotic treatment")
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
title = "after Ab treatment" if selection == "+" else "no Ab treatment"
fig.suptitle(title)
fig.tight_layout()
fig.savefig(f"./figures/LT18_GFP_calibration{selection}Ab.png",dpi=300)
fig.savefig(f"./figures/LT18_GFP_calibration{selection}Ab.svg")
with open(f"./LT_Data_py/LT18_calibration{selection}Ab.pkl", "wb") as f:
    pickle.dump(fit_result, f)
plt.show()
