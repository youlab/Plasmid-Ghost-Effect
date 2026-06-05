"""
Keio community plasmid dynamics after invasion by the plasmid-free dominant strain
Comm87+R388 <- Keio V2 #76
------------------------------------------------
• Requires:  pandas, numpy, scipy
• Input:     LT22_invasion_plating.xlsx  (sheets “Day0”, “Day2”, …,)
• Output:    plasmid_mean.npy / plasmid_std.npy
"""

import numpy as np
import pandas as pd

def mean_se(vals):
    vals = np.asarray(vals, float)
    m = np.nanmean(vals) if np.isfinite(vals).any() else np.nan
    n = np.isfinite(vals).sum()
    if n >= 2:
        sd = np.nanstd(vals, ddof=1)
        se = sd / np.sqrt(n)
    elif n == 1:
        se = np.sqrt(max(m, 0.0))  # Poisson SE
    else:
        se = np.nan
    return m, se

# ---------------------------------------------------------------------
# 1. constants & containers
# ---------------------------------------------------------------------
days           = [0, 2, 5, 8, 10, 15]
Plasmid        = "R388"

bio_rep        = 3
plate_rep      = 3

T              = len(days)
p              = 0 # p = 0 --> R388

# pre-allocate with NaN
R388_mean = np.full((bio_rep, T), np.nan)
R388_se   = np.full_like(R388_mean, np.nan)

# ---------------------------------------------------------------------
# 2. loop over sheets (“Day0”, “Day2”, …) and fill the arrays
# ---------------------------------------------------------------------
excel_file = "./raw_data/LT22_invasion_plating.xlsx"
for t_idx, t in enumerate(days):
    # read columns D:I
    # pandas keeps NaNs as float('nan') automatically
    df  = pd.read_excel(excel_file, sheet_name=f"Day{t}", usecols="C:H", skiprows=0, nrows=24)
    data = df.to_numpy(dtype=float)          # shape (6, 24)                       # plasmid condition; 0 = R388, 1 = PCU1, 2 = R6K
    for br in range(bio_rep):          # biological replicate
        # ----- locate the correct row (0-based) ----------------
        row_idx = (bio_rep * p) + br
        row_vals = data[row_idx, :2 * plate_rep]  # first 6 cols

        # LB plates = first 3 numbers, plasmid plates = next 3
        lb_vals       = row_vals[:plate_rep]
        plasmid_vals  = row_vals[plate_rep:plate_rep * 2]

        mean_LB, se_LB = mean_se(lb_vals)
        mean_plasmid, se_pl = mean_se(plasmid_vals)

        if np.isfinite(mean_LB) and mean_LB > 0 and np.isfinite(mean_plasmid):
            mean_P = 100.0 * mean_plasmid / mean_LB
            cv_LB = se_LB / mean_LB
            cv_pl = se_pl / mean_plasmid if mean_plasmid > 0 else np.nan
            se_P = mean_P * np.sqrt(cv_LB ** 2 + cv_pl ** 2)  # per-observation SE
        else:
            mean_P, se_P = np.nan, np.nan

        R388_mean[br, t_idx] = mean_P
        R388_se[br, t_idx] = se_P

# ---------------------------------------------------------------------
# 3. save one *_mean.npy and *_std.npy per plasmid
# ---------------------------------------------------------------------

np.save("./processed_data/R388_inv_mean.npy", R388_mean)
np.save("./processed_data/R388_inv_se.npy", R388_se)
