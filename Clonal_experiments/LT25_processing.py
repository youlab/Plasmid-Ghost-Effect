import numpy as np
import pandas as pd
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
plasmids = ["pSC101","colE1","pUC"]
ev_time = [3, 20]
n_clone = 9
bio_rep = 3
days = np.arange(1, 5, 1)

T = len(days)+1
for i, plasmid in enumerate(plasmids):
     for evolved_time in ev_time:
        print(f"Processing {plasmid} evolved until day {evolved_time}...")
        skip_rows = 0 if evolved_time == 3 else 3

        Mean_P = np.full((n_clone, bio_rep, T), np.nan)
        Std_P = np.full((n_clone, bio_rep, T), np.nan)

        Mean_P[:, :, 0] = 50
        Std_P[:, :, 0] = 0

        with open(f"./LT_Data_py/LT25_Day{evolved_time} {plasmid}_calibration.pkl", "rb") as f:
            fit_result = pickle.load(f)

        # Load data from Excel
        for t in days:
            # load OD and GFP raw data
            od_df = pd.read_excel(f'./Raw_data/LT25_{plasmid}.xlsx', sheet_name=f'Day{t}', header=None, skiprows=1, nrows=6,
                                usecols="B:K")
            od_np = od_df.to_numpy()
            gfp_df = pd.read_excel(f'./Raw_data/LT25_{plasmid}.xlsx', sheet_name=f'Day{t}', header=None, skiprows=9, nrows=6,
                                usecols="B:K")
            gfp_np = gfp_df.to_numpy()

            # Blank statistics
            od_blank = od_np[0:3, -1]
            gfp_blank = gfp_np[0:3, -1]

            mean_od_blank, var_od_blank = mean_and_var_of_mean(od_blank)
            mean_gfp_blank, var_gfp_blank = mean_and_var_of_mean(gfp_blank)
            # instrument/well noise (= SD of a single well, used later)
            var_od_well = np.var(od_blank, ddof=1)  # NOT divided by 3
            var_gfp_well = np.var(gfp_blank, ddof=1)

            # reading after subtracting the blank
            # clones evolved until day 3 corresponds to row A - C, until day 20 corresponds to row D - F
            od_corr = od_np[skip_rows:skip_rows+3, 0:9] - mean_od_blank
            gfp_corr = gfp_np[skip_rows:skip_rows+3, 0:9] - mean_gfp_blank

            # variance of each well’s reading = instrument + blank-mean part
            var_od_tot = var_od_well + var_od_blank
            var_gfp_tot = var_gfp_well + var_gfp_blank

            # ratio GFP/OD and its variance for every well
            y = gfp_corr / od_corr
            var_y = var_gfp_tot / od_corr ** 2 + (gfp_corr ** 2) * var_od_tot / od_corr ** 4

            for c in range(n_clone):
                p, std_p = calculate_abundance(y[:, c], var_y[:, c], fit_result[c])  # estimated plasmid abundance
                Mean_P[c,:, t] = p
                Std_P[c, :, t] = std_p

        np.save(f"./LT_data_py/LT_25_{plasmid}_day{evolved_time}_mean.npy",Mean_P)
        np.save(f"./LT_data_py/LT_25_{plasmid}_day{evolved_time}_std.npy", Std_P)
