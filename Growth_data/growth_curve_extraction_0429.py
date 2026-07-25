import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def _to_idx(well_name: list[str]) -> list[int]:
    return [well_list.index(w) for w in well_name]


colors = ["k", "#23A249",]

sample = "20260429_pUC_Day20_burden"
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

blank_wells = ["E11", "F11", "G11"]
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
    ax_gfp.set_yticks([1, 1e2, 1e3])
    ax_gfp.set_ylim([1, 4e3])
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
# fig.savefig(f"./figures/{sample}_all.png", dpi=300)

# MG1655 + pUC at Day 20
for br in range(3):
    target_wells = []
    labels = []
    for colony in range(3):
        idx = colony * 3 + br + 2
        target_wells+=[f"{r}{idx}" for r in ["B", "C", "D"]]
        labels += [f"c{colony+1}_{i}" for i in range(1,4,1)]
    br_idx = _to_idx(target_wells)
    br_OD = df_OD[br_idx, :] - blank_OD
    # Save blanked growth curve to xlsx
    OD_dict = {"time": time}
    for i, label in enumerate(labels):
        OD_dict[label] = br_OD[i, :]
    df_timeseries_od = pd.DataFrame(OD_dict)
    df_timeseries_od.to_excel(f"./processed_data/Day20/pUC+BioRep{br+1}_OD.xlsx", index=False)

# MG1655 (pUC cured) at Day 20
for br in range(3):
    target_wells = []
    labels = []
    for colony in range(3):
        idx = colony * 3 + br + 2
        target_wells+=[f"{r}{idx}" for r in ["E", "F", "G"]]
        labels += [f"c{colony+1}_{i}" for i in range(1,4,1)]
    br_idx = _to_idx(target_wells)
    br_OD = df_OD[br_idx, :] - blank_OD
    # Save blanked growth curve to xlsx
    OD_dict = {"time": time}
    for i, label in enumerate(labels):
        OD_dict[label] = br_OD[i, :]
    df_timeseries_od = pd.DataFrame(OD_dict)
    df_timeseries_od.to_excel(f"./processed_data/Day20/pUC-BioRep{br+1}_OD.xlsx", index=False)

plt.show()