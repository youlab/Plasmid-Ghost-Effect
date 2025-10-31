"""
single-strain plasmid dynamics
------------------------------------------------
• Requires:  pandas, numpy, scipy
• Input:     GE_LT10_Plating.xlsx  (sheets “Day0”, “Day3”, …, “Day25”)
• Output:    PCU1_mean.npy / PCU1_std.npy / R6K_mean.npy / R6K_std.npy
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
days           = [0, 3, 6, 7, 9, 15, 20, 25]

bio_rep        = 3
plate_rep      = 3
inits_condits  = 5
Plasmids     = ["PCU1", "R6K"]

T              = len(days)
L              = len(Plasmids)

# pre-allocate with NaN

Indiv_mean_P = np.full((L, inits_condits, bio_rep, T), np.nan)
Indiv_se_P  = np.full_like(Indiv_mean_P, np.nan)

# ---------------------------------------------------------------------
# 2. loop over sheets (“Day0”, “Day3”, …) and fill the arrays
# ---------------------------------------------------------------------
excel_file = "./Raw_data/GE_LT10_Plating.xlsx"
for t_idx, t in enumerate(days):
    # read columns D:I (0-based 3:9 exclusive) → exactly 6 columns
    # pandas keeps NaNs as float('nan') automatically
    df  = pd.read_excel(excel_file, sheet_name=f"Day{t}", usecols="D:I", skiprows=0, nrows=30)
    data = df.to_numpy(dtype=float)

    for l in range(L):                         # plasmid condition
        for ic in range(inits_condits):        # initial‐condition index
            for br in range(bio_rep):          # biological replicate
                # ----- locate the correct row (0-based) ----------------
                row_idx = (bio_rep * inits_condits * l) + (bio_rep * ic) + br
                row_vals = data[row_idx, :2 * plate_rep]              # first 6 cols

                # LB plates = first 3 numbers, plasmid plates = next 3
                lb_vals       = row_vals[:plate_rep]
                plasmid_vals  = row_vals[plate_rep:plate_rep * 2]

                mean_LB, se_LB = mean_se(lb_vals)
                mean_plasmid, se_pl = mean_se(plasmid_vals)

                if np.isfinite(mean_LB) and mean_LB > 0 and np.isfinite(mean_plasmid):
                    mean_P = 100.0 * mean_plasmid / mean_LB
                    cv_LB = se_LB / mean_LB
                    cv_pl = se_pl / mean_plasmid if mean_plasmid > 0 else 0
                    se_P = mean_P * np.sqrt(cv_LB ** 2 + cv_pl ** 2)  # per-observation SE
                else:
                    mean_P, se_P = np.nan, np.nan

                Indiv_mean_P[l, ic, br, t_idx] = mean_P
                Indiv_se_P[l,  ic, br, t_idx] = se_P

# ---------------------------------------------------------------------
# 3. save one *_mean.npy and *_std.npy per plasmid
# ---------------------------------------------------------------------
for l, plasmid in enumerate(Plasmids):
    abundance = Indiv_mean_P[l, :, :, :]         # shape (inits_condits, bio_rep, T)
    se       = Indiv_se_P[l, :, :, :]

    np.save(f"./LT_data_py/{plasmid}_mean.npy",abundance)
    np.save(f"./LT_data_py/{plasmid}_se.npy", se)
