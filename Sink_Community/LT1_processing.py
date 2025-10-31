import numpy as np
import pandas as pd

def mean_and_var_of_mean(arr):
    """Return (mean, variance_of_the_mean) for a 1-D vector."""
    n   = arr.size
    var = np.var(arr, ddof=1) / n      # divide by n → variance *of the mean*
    return arr.mean(), var

# GFP/OD level in E. coli
# ------------------------------------------------------------
# 1.  Load raw plate data
# ------------------------------------------------------------
od_df  = pd.read_excel('./raw_data/Synk_LT1.xlsx',
                       sheet_name='GFP_baseline',
                       header=None, skiprows=2,  nrows=4, usecols="B:G")
gfp_df = pd.read_excel('./raw_data/Synk_LT1.xlsx',
                       sheet_name='GFP_baseline',
                       header=None, skiprows=9,  nrows=4, usecols="B:G")
od_np  = od_df.to_numpy()
gfp_np = gfp_df.to_numpy()

# ------------------------------------------------------------
# 2.  Blank statistics (3 blank wells)
# ------------------------------------------------------------
od_blank   = od_np [-1, 0:3]
gfp_blank  = gfp_np[-1, 0:3]

mean_od_blank,  var_od_blank  = mean_and_var_of_mean(od_blank)
mean_gfp_blank, var_gfp_blank = mean_and_var_of_mean(gfp_blank)

# instrument/well noise (= SD of a single well, used later)
var_od_well  = np.var(od_blank,  ddof=1)          # NOT divided by 3
var_gfp_well = np.var(gfp_blank, ddof=1)

# ------------------------------------------------------------
# 3.  GFP/OD for pure E. coli culture (with Kan / Spect)
# ------------------------------------------------------------
od_corr  = od_np [0:3, :] - mean_od_blank
gfp_corr = gfp_np[0:3, :] - mean_gfp_blank

# variance of each well’s corrected reading = instrument + blank-mean part
var_od_tot  = var_od_well  + var_od_blank
var_gfp_tot = var_gfp_well + var_gfp_blank

# ratio GFP/OD and its SD for every well
y       = gfp_corr / od_corr
std_y   = np.sqrt(var_gfp_tot / od_corr**2 +
                  (gfp_corr**2) * var_od_tot / od_corr**4)

# Now, calculate the mean and SE of GFP/OD for each plasmid
plasmids = ["pSC101","colE1","pUC"]

for i, plasmid in enumerate(plasmids):
    baselines = np.zeros((2,2))
    yp = y[i,:]
    mean_yp_sponge = np.mean(yp[0:3])
    se_yp_sponge = np.std(yp[0:3], ddof=1) / np.sqrt(3)

    mean_yp_liquid = np.mean(yp[3:6])
    se_yp_liquid = np.std(yp[3:6], ddof=1) / np.sqrt(3)

    baselines[0,0] = mean_yp_sponge
    baselines[0,1] = se_yp_sponge

    baselines[1,0] = mean_yp_liquid
    baselines[1,1] = se_yp_liquid

    np.save(f"./processed_data/Ecoli+{plasmid}.npy",baselines)

# Initialization
days = np.arange(0, 11)
bio_rep = 3
condition = 4
plasmids = ["pSC101", "colE1", "pUC"]
T = len(days)
L = len(plasmids)

Mean_P = np.full((L, condition, bio_rep, T), np.nan)
Std_P = np.full((L, condition, bio_rep, T), np.nan)

# # fill Day 0
# Mean_P[:,:,:,0]=25
# Std_P[:,:,:,0]=0

# Load data from Excel
for t in days:
    # load OD and GFP raw data
    od_df = pd.read_excel('./raw_data/Synk_LT1.xlsx',
                          sheet_name=f'Day{t}',
                          header=None, skiprows=3, nrows=4, usecols="B:M")
    gfp_df = pd.read_excel('./raw_data/Synk_LT1.xlsx',
                           sheet_name=f'Day{t}',
                           header=None, skiprows=11, nrows=4, usecols="B:M")
    od_np = od_df.to_numpy()
    gfp_np = gfp_df.to_numpy()

    # Blank statistics
    od_blank = od_np[-1, 0:3]
    gfp_blank = gfp_np[-1, 0:3]

    mean_od_blank, var_od_blank = mean_and_var_of_mean(od_blank)
    mean_gfp_blank, var_gfp_blank = mean_and_var_of_mean(gfp_blank)
    # instrument/well noise (= SD of a single well, used later)
    var_od_well = np.var(od_blank, ddof=1)  # NOT divided by 3
    var_gfp_well = np.var(gfp_blank, ddof=1)

    # reading after subtracting the blank
    od_corr = od_np[0:3, :] - mean_od_blank
    gfp_corr = gfp_np[0:3, :] - mean_gfp_blank

    # variance of each well’s reading = instrument + blank-mean part
    var_od_tot = var_od_well + var_od_blank
    var_gfp_tot = var_gfp_well + var_gfp_blank

    # ratio GFP/OD and its variance for every well
    y = gfp_corr / od_corr
    var_y = var_gfp_tot / od_corr ** 2 + (gfp_corr ** 2) * var_od_tot / od_corr ** 4

    for l, plasmid in enumerate(plasmids):
        for c in range(condition):
            start_col = c*bio_rep
            end_col = (c+1)*bio_rep
            Mean_P[l,c,:,t] = y[l,start_col:end_col]
            Std_P[l,c,:,t] = np.sqrt(var_y[l,start_col:end_col])
# ---------------------------------------------------------------------
# save one *_mean.npy and *_std.npy per plasmid
# ---------------------------------------------------------------------
for l, plasmid in enumerate(plasmids):
    abundance = Mean_P[l, :, :, :]         # shape (conditions, bio_rep, T)
    std       = Std_P[l, :, :, :]

    np.save(f"./processed_data/{plasmid}_mean.npy",abundance)
    np.save(f"./processed_data/{plasmid}_std.npy", std)