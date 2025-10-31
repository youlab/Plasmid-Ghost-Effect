import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

def mean_and_var_of_mean(arr):
    """Return (mean, variance_of_the_mean) for a 1-D vector."""
    n   = arr.size
    var = np.var(arr, ddof=1) / n      # divide by n → variance *of the mean*
    return arr.mean(), var

def calculate_abundance(y,var_y,fit_result):
    # m is the slope, and b is the intercept
    m,b,pcov = fit_result
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
days = np.arange(1, 21)
bio_rep = 3
Ab_condits = 8
plasmids = ["pSC101", "colE1", "pUC"]
T = len(days)+1
L = len(plasmids)

Mean_P = np.full((L, Ab_condits, bio_rep, T), np.nan)
Std_P = np.full((L, Ab_condits, bio_rep, T), np.nan)

Mean_P[:,:,:,0] = 50
Std_P[:,:,:,0] = 0

# Load data from Excel
for t in days:
    with open(f"./LT_Data_py/LT18_calibration+Ab.pkl", "rb") as f:
        fit_result = pickle.load(f)
    # load OD and GFP raw data
    od_df = pd.read_excel('./Raw_data/GE_LT18_GFP.xlsx', sheet_name=f'Day{t}', header=None, skiprows=2, nrows=8,
                          usecols="B:M")
    od_np = od_df.to_numpy()
    gfp_df = pd.read_excel('./Raw_data/GE_LT18_GFP.xlsx', sheet_name=f'Day{t}', header=None, skiprows=13, nrows=8,
                          usecols="B:M")
    gfp_np = gfp_df.to_numpy()

    # Blank statistics
    od_blank = od_np[4, 9:12]
    gfp_blank = gfp_np[4, 9:12]

    mean_od_blank, var_od_blank = mean_and_var_of_mean(od_blank)
    mean_gfp_blank, var_gfp_blank = mean_and_var_of_mean(gfp_blank)
    # instrument/well noise (= SD of a single well, used later)
    var_od_well = np.var(od_blank, ddof=1)  # NOT divided by 3
    var_gfp_well = np.var(gfp_blank, ddof=1)

    # reading after subtracting the blank
    od_corr = od_np[:, 0:9] - mean_od_blank
    gfp_corr = gfp_np[:, 0:9] - mean_gfp_blank

    # variance of each well’s reading = instrument + blank-mean part
    var_od_tot = var_od_well + var_od_blank
    var_gfp_tot = var_gfp_well + var_gfp_blank

    # ratio GFP/OD and its variance for every well
    y = gfp_corr / od_corr
    var_y = var_gfp_tot / od_corr ** 2 + (gfp_corr ** 2) * var_od_tot / od_corr ** 4

    for l, plasmid in enumerate(plasmids):
        start_col = l * bio_rep
        end_col = (l + 1) * bio_rep
        p, std_p = calculate_abundance(y[:,start_col:end_col], var_y[:,start_col:end_col], fit_result[plasmid])# estimated plasmid abundance
        Mean_P[l, :, :, t] = p
        Std_P[l, :, :, t] = std_p
# ---------------------------------------------------------------------
# save one *_mean.npy and *_std.npy per plasmid
# ---------------------------------------------------------------------
for l, plasmid in enumerate(plasmids):
    abundance = Mean_P[l, :, :, :]         # shape (inits_condits, bio_rep, T)
    std       = Std_P[l, :, :, :]

    np.save(f"./LT_data_py/{plasmid}_dose_response_mean.npy",abundance)
    np.save(f"./LT_data_py/{plasmid}_dose_response_std.npy", std)
