import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def _to_idx(well_name: list[str]) -> list[int]:
    return [well_list.index(w) for w in well_name]


colors = ["k", "#23A249", "#B53030", "#74B3EB", "#FF9000"]

sample = "20260318_GE_Ancestor_burden"
file = f"./raw_data/{sample}.xlsx"
time = pd.read_excel(file, sheet_name="Sheet2", header=None, usecols="B:KC", skiprows=47, nrows=1)
time = time.to_numpy(dtype=float)[0, :] / 3600 # convert to hours

ncol_plate = 12 # rows in the 96-well plate
nrow_plate = 8 # columns in the 96-well plate

header_wells = pd.read_excel(file, sheet_name="Sheet2",header=None, usecols="A", skiprows=49, nrows=96)
well_list = header_wells.iloc[:, 0].astype(str).tolist()

# load data
df_OD = pd.read_excel(file, sheet_name="Sheet2", header=None, usecols="B:KC", skiprows=49, nrows=96).to_numpy(dtype=float)
df_GFP = pd.read_excel(file, sheet_name="Sheet2", header=None, usecols="B:KC", skiprows=150, nrows=96).to_numpy(dtype=float)

blank_wells = ["B11", "C11", "D11"]
blank_idx = _to_idx(blank_wells)

blank_OD = np.mean(df_OD[blank_idx, :], axis=0)
blank_GFP = np.mean(df_GFP[blank_idx, :], axis=0)

fig, axes = plt.subplots(nrow_plate, ncol_plate, figsize=(0.8*ncol_plate, 0.8*nrow_plate))
axes = axes.flat
for c in range(ncol_plate*nrow_plate):
    ax = axes[c]
    ax_gfp = ax.twinx()  # shares x with ax; separate right y-axis for GFP
    OD_raw = df_OD[c, :]
    GFP_raw = df_GFP[c, :]
    ax.plot(time, OD_raw, linewidth=1.2, c=colors[0])
    ax_gfp.plot(time, GFP_raw, linewidth=1.2, c=colors[1])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax_gfp.set_yscale("log")
    ax_gfp.set_yticks([1, 1e2, 1e4])
    ax_gfp.set_ylim([1, 1e5])
    ax.set_xticks([0, 12, 24])
    if c == ncol_plate*(nrow_plate - 1):
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("OD600")
        ax_gfp.set_yticklabels([])
    elif c == ncol_plate*nrow_plate - 1:
        ax_gfp.set_ylabel("GFP", color=colors[1])
        ax.set_yticklabels([])
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax_gfp.set_yticklabels([])
fig.suptitle(sample)
fig.subplots_adjust(hspace=0.05, wspace=0.05,left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.savefig(f"./figures/{sample}_all.png", dpi=300)

# pSC101: LB + Kan
target_wells = ["B2","B3","B4","C2","C3","C4","D2","D3","D4"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
pSC101_idx = _to_idx(target_wells)
pSC101_OD = df_OD[pSC101_idx, :] - blank_OD
pSC101_GFP = df_GFP[pSC101_idx, :] - blank_GFP

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
GFP_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = pSC101_OD[i, :]
    GFP_dict[label] = pSC101_GFP[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/pSC101_OD.xlsx", index=False)

# colE1: LB + Kan
target_wells = ["B5","B6","B7","C5","C6","C7","D5","D6","D7"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
colE1_idx = _to_idx(target_wells)
colE1_OD = df_OD[colE1_idx, :] - blank_OD
colE1_GFP = df_GFP[colE1_idx, :] - blank_GFP

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
GFP_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = colE1_OD[i, :]
    GFP_dict[label] = colE1_GFP[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/colE1_OD.xlsx", index=False)

# pUC: LB + Spect
target_wells = ["B8","B9","B10","C8","C9","C10","D8","D9","D10"]
labels = ["c1_1","c2_1","c3_1","c1_2","c2_2","c3_2","c1_3","c2_3","c3_3"]
pUC_idx = _to_idx(target_wells)
pUC_OD = df_OD[pUC_idx, :] - blank_OD
pUC_GFP = df_GFP[pUC_idx, :] - blank_GFP

# Save blanked growth curve to xlsx
OD_dict = {"time": time}
GFP_dict = {"time": time}
for i, label in enumerate(labels):
    OD_dict[label] = pUC_OD[i, :]
    GFP_dict[label] = pUC_GFP[i, :]

df_timeseries_od = pd.DataFrame(OD_dict)
df_timeseries_od.to_excel(f"./processed_data/ancestor/pUC_OD.xlsx", index=False)
plt.show()