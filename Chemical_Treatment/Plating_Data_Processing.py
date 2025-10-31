"""
pSC101 dynamics with chemical treatments
------------------------------------------------
• Requires:  pandas, numpy, scipy
• Input:     Chemical_LT14.xlsx  (sheets “Day0”, “Day1”, …, “Day10”)
• Output:    Rif_mean.npy / Rif_std.npy / Rif+Prm_mean.npy / Rif+Prm_std.npy / Rif+Pmt_mean.npy / Rif+Pmt_std.npy / Rif+Phe_mean.npy / Rif+Phe_std.npy / NoChem_mean.npy / NoChem_std.npy
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
days           = [0, 1, 2, 3, 4, 5, 10]

bio_rep        = 3
plate_rep      = 3
treatment_durations  = 4 # 1-day, 2-day, 3-day, and continuous
Chemicals     = ["Rif", "Rif+Prm", "Rif+Mpt", "Rif+Phe"]

T              = len(days)
L              = len(Chemicals)

# pre-allocate with NaN

Indiv_mean_P = np.full((L, treatment_durations, bio_rep, T), np.nan)
Indiv_se_P  = np.full_like(Indiv_mean_P, np.nan)

Chem_free_mean = np.full((bio_rep,T),np.nan)
Chem_free_se = np.full_like(Chem_free_mean, np.nan)

# ---------------------------------------------------------------------
# 2. loop over sheets (“Day0”, “Day1”, …) and fill the arrays
# ---------------------------------------------------------------------
excel_file = "./Chemical_LT14.xlsx"
for t_idx, t in enumerate(days):
    # read columns D:I (0-based 3:9 exclusive) → exactly 6 columns
    # pandas keeps NaNs as float('nan') automatically
    df  = pd.read_excel(excel_file, sheet_name=f"Day{t}", usecols="D:I", skiprows=0, nrows=51)
    data = df.to_numpy(dtype=float)          # shape (51, 6)

    for l in range(L):                         # chemical condition
        for tc in range(treatment_durations):        # treatment duration index (1 / 2 / 3 -day treatment)
            for br in range(bio_rep):          # biological replicate
                # ----- locate the correct row (0-based) ----------------
                row_idx = (bio_rep * L * tc) + (bio_rep * l) + br
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

                Indiv_mean_P[l, tc, br, t_idx] = mean_P
                Indiv_se_P[l,  tc, br, t_idx] = se_P
    # for chemical-free condition
    for br in range(3):
        cf_idx = bio_rep * L * treatment_durations + br
        row_vals = data[cf_idx, :2 * plate_rep]
        # LB plates = first 3 numbers, plasmid plates = next 3
        lb_vals = row_vals[:plate_rep]
        plasmid_vals = row_vals[plate_rep:plate_rep * 2]

        mean_LB, se_LB = mean_se(lb_vals)
        mean_plasmid, se_pl = mean_se(plasmid_vals)

        if np.isfinite(mean_LB) and mean_LB > 0 and np.isfinite(mean_plasmid):
            mean_P = 100.0 * mean_plasmid / mean_LB
            cv_LB = se_LB / mean_LB
            cv_pl = se_pl / mean_plasmid if mean_plasmid > 0 else 0
            se_P = mean_P * np.sqrt(cv_LB ** 2 + cv_pl ** 2)  # per-observation SE
        else:
            mean_P, se_P = np.nan, np.nan
        Chem_free_mean[br, t_idx] = mean_P
        Chem_free_se[br, t_idx] = se_P
# ---------------------------------------------------------------------
# 3. save one *_mean.npy and *_std.npy per plasmid
# ---------------------------------------------------------------------
for l, chem in enumerate(Chemicals):
    abundance = Indiv_mean_P[l, :, :, :]         # shape (treatment_durations, bio_rep, T)
    se       = Indiv_se_P[l, :, :, :]
    se[:,:,0] = np.nan

    np.save(f"./processed_data/{chem}_mean.npy",abundance)
    np.save(f"./processed_data/{chem}_se.npy", se)
Chem_free_se[:,0] = np.nan
np.save(f"./processed_data/NoChem_mean.npy", Chem_free_mean)
np.save(f"./processed_data/NoChem_se.npy", Chem_free_se)
