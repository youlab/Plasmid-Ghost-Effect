import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def _to_idx(well_name: list[str]) -> list[int]:
    return [well_list.index(w) for w in well_name]


colors = ["k", "#23A249", "#B53030", "#74B3EB", "#FF9000"]

sample = "20260323_GE_DA_Ancestor_burden"
file = f"./raw_data/{sample}.xlsx"
time = pd.read_excel(file, sheet_name="Sheet2", header=None, usecols="B:KC", skiprows=47, nrows=1)
time = time.to_numpy(dtype=float)[0, :] / 3600 # convert to hours

ncol_plate = 12 # rows in the 96-well plate
nrow_plate = 8 # columns in the 96-well plate

header_wells = pd.read_excel(file, sheet_name="Sheet2",header=None, usecols="A", skiprows=49, nrows=96)
well_list = header_wells.iloc[:, 0].astype(str).tolist()

# load data
df_OD = pd.read_excel(file, sheet_name="Sheet2", header=None, usecols="B:KC", skiprows=49, nrows=96).to_numpy(dtype=float)

blank_wells = ["E11", "F11", "G11"]
blank_idx = _to_idx(blank_wells)

blank_OD = np.mean(df_OD[blank_idx, :], axis=0)

fig, axes = plt.subplots(nrow_plate, ncol_plate, figsize=(0.8*ncol_plate, 0.8*nrow_plate))
axes = axes.flat
for c in range(ncol_plate*nrow_plate):
    ax = axes[c]
    OD_raw = df_OD[c, :]
    ax.plot(time, OD_raw, linewidth=1.2, c=colors[0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xticks([0, 12, 24])
    if c == ncol_plate*(nrow_plate - 1):
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("OD600")
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
fig.suptitle(sample)
fig.subplots_adjust(hspace=0.05, wspace=0.05,left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.savefig(f"./figures/{sample}_all.png", dpi=300)

# DA28102: LB + Cm
target_wells = ["B2","B3","B4","C2","C3","C4","D2","D3","D4"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
DA_idx = _to_idx(target_wells)
DA_OD = df_OD[DA_idx, :] - blank_OD

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = DA_OD[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/DA28102_OD.xlsx", index=False)

# MG1655: LB
target_wells = ["E2","E3","E4","F2","F3","F4","G2","G3","G4"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
MG_idx = _to_idx(target_wells)
MG_OD = df_OD[MG_idx, :] - blank_OD

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
GFP_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = MG_OD[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/MG1655_OD.xlsx", index=False)

# pCU1: LB + Cm + Carb
target_wells = ["E5","E6","E7","F5","F6","F7","G5","G6","G7"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
pCU1_idx = _to_idx(target_wells)
pCU1_OD = df_OD[pCU1_idx, :] - blank_OD

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = pCU1_OD[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/pCU1_OD.xlsx", index=False)

# R6K: LB + Cm + Strp
target_wells = ["E8","E9","E10","F8","F9","F10","G8","G9","G10"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
R6K_idx = _to_idx(target_wells)
R6K_OD  = df_OD[R6K_idx, :] - blank_OD

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = R6K_OD[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/R6K_OD.xlsx", index=False)
plt.show()