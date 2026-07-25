# processing qPCR data for MG1655 + pSC101
# this file includes loading calibration data, and generate estimates for individual samples
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t
from sklearn.metrics import r2_score

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

# Linear model
def linear_model(x, m, b):
    return m * x + b

def linear_fit(x, y):
    """Perform linear regression and return slope, intercept, and covariance."""
    y_nan = np.isnan(y)
    x = x[~y_nan]
    y = y[~y_nan]
    popt, pcov = curve_fit(linear_model, x, y)
    slope, intercept = popt
    slope_se, intercept_se = np.sqrt(np.diag(pcov))
    # Compute t-statistic and p-value
    dof = len(x) - 2
    t_stat = slope / slope_se
    p_value = 2 * (1 - t.cdf(np.abs(t_stat), df=dof))
    # residuals and total sums of squares
    y_pred = linear_model(x, slope, intercept)

    R2 = r2_score(y, y_pred)

    # --- 95 % confidence band -----------------------------------------------
    t_crit = t.ppf(0.975, dof)  # two-tailed 95 %

    x_band = np.linspace(np.min(x), np.max(x), 300)  # dense grid for smooth band
    y_band = linear_model(x_band, slope, intercept)

    # Jacobian of y wrt [slope, intercept] at every x -> [[x, 1], ...]
    J = np.vstack((x_band, np.ones_like(x_band))).T  # shape (N, 2)

    # point-wise variance:  diag( J · pcov · Jᵀ )
    sigma2_y = np.einsum('ij,jk,ik->i', J, pcov, J)
    ci_y = t_crit * np.sqrt(sigma2_y)  # half-width of CI
    return popt, pcov, p_value, R2, x_band, y_band, ci_y

# Calibration data
file = "./Raw_data/pSC101_colE1_calibration.xlsx"

df = pd.read_excel(file, sheet_name="Results")

df_fam = df[df["Reporter"] == "FAM"].reset_index(drop=True)
df_vic = df[df["Reporter"] == "VIC"].reset_index(drop=True)

# MG1655 + pSC101
log10_dilution = np.arange(0, -8, -1)  # dilution: 10^0, 10^-1, ..., 10^-7
log2_dilution = np.log2(10.0 ** log10_dilution)  # convert to log2 scale
pSC_fam = np.zeros((8, 6)) * np.nan
pSC_vic = np.zeros((8, 6)) * np.nan

for i, row in enumerate("ABCDEFGH"):
    for j, col in enumerate(range(1, 7)):
        well = f"{row}{col}"
        fam_Ct = df_fam.loc[df_fam["Well Position"] == well, "CT"].values[0]
        vic_Ct = df_vic.loc[df_vic["Well Position"] == well, "CT"].values[0]
        if fam_Ct == "Undetermined":
            fam_Ct = np.nan
        if vic_Ct == "Undetermined":
            vic_Ct = np.nan
        pSC_fam[i, j] = fam_Ct
        pSC_vic[i, j] = vic_Ct

# plot calibration curve
fig,ax = plt.subplots(figsize=(3, 3))
for i in range(6):
    ax.scatter(log2_dilution, pSC_fam[:,i], color=colors[0], marker="o", edgecolor=colors[0], facecolor="None", s=50, lw=1.5)
    ax.scatter(log2_dilution, pSC_vic[:,i], color=colors[1], marker="o", edgecolor=colors[1], facecolor="None", s=50, lw=1.5)   

# Weighted fit 
x = np.tile(log2_dilution, 6)  # repeat dilution for each replicate
popt_fam, pcov_fam, p_fam, R2_fam, x_band_fam, y_band_fam, ci_y_fam = linear_fit(x, pSC_fam.T.flatten())
popt_vic, pcov_vic, p_vic, R2_vic, x_band_vic, y_band_vic, ci_y_vic = linear_fit(x, pSC_vic.T.flatten())
fam_slope_se = np.sqrt(pcov_fam[0, 0])
vic_slope_se = np.sqrt(pcov_vic[0, 0])

slope_summary = np.array([
    [popt_fam[0], fam_slope_se],
    [popt_vic[0], vic_slope_se],
])
np.savetxt(
    "./processed_data/MG1655-pSC101-calibration.txt",
    slope_summary,
    fmt="%.6f",
    delimiter="\t",
)

# plot fit line and its confidence envelope
ax.plot(x_band_fam, y_band_fam, c=colors[0], lw=1)
ax.fill_between(x_band_fam, y_band_fam - ci_y_fam, y_band_fam + ci_y_fam,
                color='steelblue', alpha=0.50, linewidth=0)
ax.plot(x_band_vic, y_band_vic, c=colors[1], lw=1)
ax.fill_between(x_band_vic, y_band_vic - ci_y_vic, y_band_vic + ci_y_vic,
                color='steelblue', alpha=0.50, linewidth=0)
# Format regression equation
txt = f"m = {popt_fam[0]:.3f} ± {fam_slope_se:.3f}\nR$^2$ = {R2_fam:.2f}"
ax.text(x=0.9, y=0.78, s=txt, transform=ax.transAxes, fontsize=12, ha='right', c= colors[0])
txt = f"m = {popt_vic[0]:.3f} ± {vic_slope_se:.3f}\nR$^2$ = {R2_vic:.2f}"
ax.text(x=0.1, y=0.1, s=txt, transform=ax.transAxes, fontsize=12, c = colors[1])
ax.set_xlabel("log$_{2}$dilution")
ax.set_ylabel("Ct")
ax.set_title("MG1655+pSC101")
fig.tight_layout()
# fig.savefig("./figures/MG1655-pSC101-calibration.png", dpi=300)
# fig.savefig("./figures/MG1655-pSC101-calibration.svg", dpi=300)

# measurement data
file = "./Raw_data/measurements.xlsx"

df = pd.read_excel(file, sheet_name="Results")

df_fam = df[df["Reporter"] == "FAM"].reset_index(drop=True)
df_vic = df[df["Reporter"] == "VIC"].reset_index(drop=True)

# MG1655 + pSC101
pSC_fam = np.zeros((3, 3)) * np.nan # biological replicates * technical replicates
pSC_vic = np.zeros((3, 3)) * np.nan

for i, row in enumerate("ABC"):
    for j, col in enumerate(range(7, 10)):
        well = f"{row}{col}"
        fam_Ct = df_fam.loc[df_fam["Well Position"] == well, "CT"].values[0]
        vic_Ct = df_vic.loc[df_vic["Well Position"] == well, "CT"].values[0]
        if fam_Ct == "Undetermined":
            fam_Ct = np.nan
        if vic_Ct == "Undetermined":
            vic_Ct = np.nan
        pSC_fam[i, j] = fam_Ct
        pSC_vic[i, j] = vic_Ct

CN = 2**(pSC_vic/popt_vic[0]) / 2**(pSC_fam/popt_fam[0]) / 2 # plasmid pSC101 is a concatemer of 2 copies, so divide by 2 to get copy number per plasmid molecule
np.savetxt("./processed_data/pSC101_CN.txt",CN)