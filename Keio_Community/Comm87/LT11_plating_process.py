"""
Keio V2 community plasmid dynamics
------------------------------------------------
• Requires:  pandas, numpy, scipy
• Input:     GE_LT11_Keio.xlsx  (sheets “Day0”, “Day3”, …, “Day25”)
• Output:    R388_mean.npy / R388_std.npy
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
days           = [0, 2, 3, 7, 10, 15, 20, 25, 30, 36]
plasmid        = "R388"

bio_rep        = 3
plate_rep      = 3
Ab_condits     = 2

T              = len(days)

# pre-allocate with NaN

Indiv_mean_P = np.full((Ab_condits, bio_rep, T), np.nan)
Indiv_se_P  = np.full_like(Indiv_mean_P, np.nan)

# ---------------------------------------------------------------------
# 2. loop over sheets (“Day0”, “Day2”, …) and fill the arrays
# ---------------------------------------------------------------------
excel_file = "./raw_data/GE_LT11_Keio.xlsx"
for t_idx, t in enumerate(days):
    # read columns C:H
    # pandas keeps NaNs as float('nan') automatically
    df  = pd.read_excel(excel_file, sheet_name=f"Day{t}", usecols="C:H", skiprows=0, nrows=6)
    data = df.to_numpy(dtype=float)          # shape (6, 6)
    for ab in range(Ab_condits):        # Ab pulsed vs no pulse
        for br in range(bio_rep):          # biological replicate
            # ----- locate the correct row (0-based) ----------------
            row_idx = (bio_rep * ab) + br
            row_vals = data[row_idx, :2 * plate_rep]              # first 6 cols

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

            Indiv_mean_P[ab, br, t_idx] = mean_P
            Indiv_se_P[ab, br, t_idx] = se_P

# ---------------------------------------------------------------------
# 3. save one R388_mean.npy and R388_std.npy
# ---------------------------------------------------------------------

np.save(f"./processed_data/{plasmid}_mean.npy",Indiv_mean_P)
np.save(f"./processed_data/{plasmid}_se.npy", Indiv_se_P)
