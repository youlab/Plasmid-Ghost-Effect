import pandas as pd
import re
import os
import numpy as np

# ===== Settings =====
input_file = "./raw_data/16S_composition.xlsx"   # change to your file path
output_dir = "./processed_data/"                # folder to save outputs
taxonomy_column = "Amplicon Sequencing Variant (ASV)"  # adjust if needed
rank_prefix = "g__"                  # collapse to genus
cutoff = 0.01                        # 1% cutoff (use 0.001 for 0.1%)

# ===== Functions =====
def extract_rank(tax_str: str, rank_prefix: str = "g__") -> str:
    """Extracts the requested rank (e.g., genus 'g__') from taxonomy string."""
    if not isinstance(tax_str, str):
        return "Unclassified"
    parts = [p.strip() for p in tax_str.split(";") if p.strip()]
    for p in parts:
        if p.startswith(rank_prefix):
            name = p.replace(rank_prefix, "", 1)
            return name if name else "Unclassified"
    labeled = [re.sub(r"^[a-z]__", "", p) for p in parts if p]
    return labeled[-1] if labeled else "Unclassified"

# ===== Load data =====
df = pd.read_excel(input_file, sheet_name="16S_composition")

# Identify sample columns = numeric columns other than taxonomy
sample_cols = [c for c in df.columns if c != taxonomy_column and pd.api.types.is_numeric_dtype(df[c])]

# Extract genus and collapse (this stays pristine)
df["Genus"] = df[taxonomy_column].apply(lambda s: extract_rank(s, rank_prefix=rank_prefix))
genus_table = df.groupby("Genus")[sample_cols].sum().reset_index()

# ===== Build plot-ready table: abundant genera + "Other" =====
# 1) Determine which genera are "abundant" anywhere
abundant_mask_rows = (genus_table[sample_cols].max(axis=1) >= cutoff)
abundant = genus_table.loc[abundant_mask_rows].copy()

# 2) In a *copy* used for plotting, zero-out sub-cutoff cells for abundant genera
#    (so those values can be moved into "Other" while keeping abundant ones where they matter)
low_cells_mask = abundant[sample_cols] < cutoff
abundant_zeroed = abundant.copy()
abundant_zeroed[sample_cols] = abundant_zeroed[sample_cols].where(~low_cells_mask, 0)

# 3) Build the "Other" row: sum of *all* sub-cutoff cells across ALL genera
#    (includes: a) entirely rare genera dropped from 'abundant', and b) the <cutoff cells
#     of abundant genera that were zeroed above)
#    To get that, compute column-wise totals from the full table and subtract the kept abundant_zeroed totals.
full_totals = genus_table[sample_cols].sum(axis=0)
kept_totals = abundant_zeroed[sample_cols].sum(axis=0)
other_values = (full_totals - kept_totals).tolist()

other_row = pd.DataFrame([["Other"] + other_values], columns=["Genus"] + sample_cols)

# 4) Concatenate abundant_zeroed + Other
genus_filtered = pd.concat([abundant_zeroed, other_row], ignore_index=True)

# (Optional) If any sample has tiny FP drift, you can normalize here; not required since nothing was discarded.

# ===== Save outputs =====
os.makedirs(output_dir, exist_ok=True)
out_all = os.path.join(output_dir, "genus_level_composition.xlsx")                # pristine
out_filtered = os.path.join(output_dir, "genus_level_composition_filtered.xlsx")  # plot-ready

genus_table.to_excel(out_all, index=False)
genus_filtered.to_excel(out_filtered, index=False)

print(f"[ok] Collapsed to genus level (pristine) → {out_all}")
print(f"[ok] Plot-ready (abundant + Other, cutoff < {cutoff:.2%}) → {out_filtered}")
print(f"[info] Abundant genera kept: {abundant.shape[0]} (plus 'Other')")

# ===== Construct numpy arrays each with a particular biological replicate for plotting later =====
# 1) Define the samples between Day0 - 2 for each biological condition
groups = {
    "pSC101": {
        "NoAb": [
            ["A1","A4","B1"],
            ["A1","A5","B2"],
            ["A1","A6","B3"]
        ],
        "Ab": [
            ["A1","A4","C1"],
            ["A1","A5","C2"],
            ["A1","A6","C3"]
        ]
    },
    "colE1": {
        "NoAb": [
            ["A2","A7","B4"],
            ["A2","A8","B5"],
            ["A2","A9","B6"],
        ],
        "Ab": [
            ["A2","A7","C4"],
            ["A2","A8","C5"],
            ["A2","A9","C6"],
        ]
    },
    "pUC": {
        "NoAb": [
            ["A3","A10","B7"],
            ["A3","A11","B8"],
            ["A3","A12","B9"]
        ],
        "Ab": [
            ["A3","A10","C7"],
            ["A3","A11","C8"],
            ["A3","A12","C9"]
        ]
    }
}


# ===== Load data =====
df = pd.read_excel(out_filtered)
print(df)

plasmids = ["pSC101","colE1","pUC"]
Ab_conditions = ["NoAb","Ab"]

for i,plasmid in enumerate(plasmids):
    for j,cond in enumerate(Ab_conditions):
        replicates = groups[plasmid][cond]
        data = []
        for rep in replicates:
            subset = df[rep].to_numpy()
            data.append(subset)
        data = np.array(data)
        np.save(f"./processed_data/{plasmid}_{cond}.npy",data)

