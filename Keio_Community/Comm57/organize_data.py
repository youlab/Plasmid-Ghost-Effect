import pandas as pd
import numpy as np
def save_data(data, dict_list, sample_name, save=True):
    """
    Splits the data into multiple DataFrames based on the provided column-to-sample dictionaries
    and renames the columns accordingly.

    Args:
        data (pd.DataFrame): Original DataFrame.
        dict_list (list): List of dictionaries mapping columns to sample names.
        sample_name (str): Prefix for the output filenames.

    Returns:
        list: A list of DataFrames, each corresponding to a dictionary.
    """
    dfs = []
    antibiotics = ["NoAb","Ab"]
    with pd.ExcelWriter("./organized_NGS/%s.xlsx"%sample_name) as writer:
        for i,column_to_sample_dicts in enumerate(dict_list):
            ab_label = antibiotics[i]
            for j, column_to_sample in enumerate(column_to_sample_dicts):

                rep = j+1
                # Filter the columns that are in the current dictionary
                sorted_columns = sorted(
                    column_to_sample.items(),
                    key=lambda x: int(x[1].split("Day ")[1])  # Extract numeric day and sort
                )
                sorted_column_names = ["Barcodes",]+[col[0] for col in sorted_columns]

                # Filter and reorder the columns in the DataFrame
                df_split = data[sorted_column_names].copy()

                # Rename the columns based on the sorted dictionary
                df_split.rename(columns=dict(sorted_columns), inplace=True)

                # Save the DataFrame to the current sheet
                if save:
                    df_split.to_excel(writer, sheet_name="%s_Rep%i"%(ab_label,rep), index=False)
    if save:
        print("%s excel file saved"%sample_name)

save = True
# Step 1: Read the CSV file
file_path = "merged_counts.csv"  # Replace with your actual file path
data = pd.read_csv(file_path)

# Display the first few rows to understand the format
print("Original Data:")
print(data.head())
# Step 2: Define multiple dictionaries to map columns to sample names
# Example mappings; modify these as needed
"""
1:500 dilution
"""
# Keio #1 + R388
# three replicates; two conditions: with or without antibiotic pulse
R388_500_no_pulse = [
    {
        "1-1": "Day 0",
        "1-3": "Day 2",
        "1-7": "Day 3",
        "1-13": "Day 7",
        "1-19": "Day 10",
        "9-1": "Day 15",
        "9-7": "Day 20",
    },
    {
        "1-1": "Day 0",
        "2-3": "Day 2",
        "1-8": "Day 3",
        "1-14": "Day 7",
        "1-20": "Day 10",
        "9-2": "Day 15",
        "9-8": "Day 20",
    },
    {
        "1-1": "Day 0",
        "3-3": "Day 2",
        "1-9": "Day 3",
        "1-15": "Day 7",
        "1-21": "Day 10",
        "9-3": "Day 15",
        "9-9": "Day 20",
    },
]
R388_500_pulse = [
    {
        "1-1": "Day 0",
        "1-3": "Day 2",
        "5-7": "Day 3",
        "5-13": "Day 7",
        "5-19": "Day 10",
        "13-1": "Day 15",
        "13-7": "Day 20",
    },
    {
        "1-1": "Day 0",
        "2-3": "Day 2",
        "5-8": "Day 3",
        "5-14": "Day 7",
        "5-20": "Day 10",
        "13-2": "Day 15",
        "13-8": "Day 20",
    },
    {
        "1-1": "Day 0",
        "3-3": "Day 2",
        "5-9": "Day 3",
        "5-15": "Day 7",
        "5-21": "Day 10",
        "13-3": "Day 15",
        "13-9": "Day 20",
    },
]
R388_dicts = [R388_500_no_pulse,R388_500_pulse]
save_data(data,R388_dicts,"R388_500_dilution",save=save)
# Keio #2 + RP4
RP4_500_no_pulse = [
    {
        "2-1": "Day 0",
        "1-4": "Day 2",
        "2-7": "Day 3",
        "2-13": "Day 7",
        "2-19": "Day 10",
        "10-1": "Day 15",
        "10-7": "Day 20",
    },
    {
        "2-1": "Day 0",
        "2-4": "Day 2",
        "2-8": "Day 3",
        "2-14": "Day 7",
        "2-20": "Day 10",
        "10-2": "Day 15",
        "10-8": "Day 20",
    },
    {
        "2-1": "Day 0",
        "3-4": "Day 2",
        "2-9": "Day 3",
        "2-15": "Day 7",
        "2-21": "Day 10",
        "10-3": "Day 15",
        "10-9": "Day 20",
    },
]
RP4_500_pulse = [
    {
        "2-1": "Day 0",
        "1-4": "Day 2",
        "6-7": "Day 3",
        "6-13": "Day 7",
        "6-19": "Day 10",
        "14-1": "Day 15",
        "14-7": "Day 20",
    },
    {
        "2-1": "Day 0",
        "2-4": "Day 2",
        "6-8": "Day 3",
        "6-14": "Day 7",
        "6-20": "Day 10",
        "14-2": "Day 15",
        "14-8": "Day 20",
    },
    {
        "2-1": "Day 0",
        "3-4": "Day 2",
        "6-9": "Day 3",
        "6-15": "Day 7",
        "6-21": "Day 10",
        "14-3": "Day 15",
        "14-9": "Day 20",
    },
]
RP4_dicts = [RP4_500_no_pulse,RP4_500_pulse]
save_data(data,RP4_dicts,"RP4_500_dilution",save=save)
# Keio #3 + pCU1
pCU1_500_no_pulse = [
    {
        "3-1": "Day 0",
        "1-5": "Day 2",
        "3-7": "Day 3",
        "3-13": "Day 7",
        "3-19": "Day 10",
        "11-1": "Day 15",
        "11-7": "Day 20",
    },
    {
        "3-1": "Day 0",
        "2-5": "Day 2",
        "3-8": "Day 3",
        "3-14": "Day 7",
        "3-20": "Day 10",
        "11-2": "Day 15",
        "11-8": "Day 20",
    },
    {
        "3-1": "Day 0",
        "3-5": "Day 2",
        "3-9": "Day 3",
        "3-15": "Day 7",
        "3-21": "Day 10",
        "11-3": "Day 15",
        "11-9": "Day 20",
    },
]
pCU1_500_pulse = [
    {
        "3-1": "Day 0",
        "1-5": "Day 2",
        "7-7": "Day 3",
        "7-13": "Day 7",
        "7-19": "Day 10",
        "15-1": "Day 15",
        "15-7": "Day 20",
    },
    {
        "3-1": "Day 0",
        "2-5": "Day 2",
        "7-8": "Day 3",
        "7-14": "Day 7",
        "7-20": "Day 10",
        "15-2": "Day 15",
        "15-8": "Day 20",
    },
    {
        "3-1": "Day 0",
        "3-5": "Day 2",
        "7-9": "Day 3",
        "7-15": "Day 7",
        "7-21": "Day 10",
        "15-3": "Day 15",
        "15-9": "Day 20",
    },
]
pCU1_dicts = [pCU1_500_no_pulse,pCU1_500_pulse]
save_data(data,pCU1_dicts,"pCU1_500_dilution",save=save)
# Keio #4 + R6K
R6K_500_no_pulse = [
    {
        "4-1": "Day 0",
        "1-6": "Day 2",
        "4-7": "Day 3",
        "4-13": "Day 7",
        "4-19": "Day 10",
        "12-1": "Day 15",
        "12-7": "Day 20",
    },
    {
        "4-1": "Day 0",
        "2-6": "Day 2",
        "4-8": "Day 3",
        "4-14": "Day 7",
        "4-20": "Day 10",
        "12-2": "Day 15",
        "12-8": "Day 20",
    },
    {
        "4-1": "Day 0",
        "3-6": "Day 2",
        "4-9": "Day 3",
        "4-15": "Day 7",
        "4-21": "Day 10",
        "12-3": "Day 15",
        "12-9": "Day 20",
    },
]
R6K_500_pulse = [
    {
        "4-1": "Day 0",
        "1-6": "Day 2",
        "8-7": "Day 3",
        "8-13": "Day 7",
        "8-19": "Day 10",
        "16-1": "Day 15",
        "16-7": "Day 20",
    },
    {
        "4-1": "Day 0",
        "2-6": "Day 2",
        "8-8": "Day 3",
        "8-14": "Day 7",
        "8-20": "Day 10",
        "16-2": "Day 15",
        "16-8": "Day 20",
    },
    {
        "4-1": "Day 0",
        "3-6": "Day 2",
        "8-9": "Day 3",
        "8-15": "Day 7",
        "8-21": "Day 10",
        "16-3": "Day 15",
        "16-9": "Day 20",
    },
]
R6K_dicts = [R6K_500_no_pulse,R6K_500_pulse]
save_data(data,R6K_dicts,"R6K_500_dilution",save=save)

# confirm the dictionaries are correct
import matplotlib.pyplot as plt
import seaborn as sns
all_dicts = R388_dicts+RP4_dicts+pCU1_dicts+R6K_dicts
Plasmids = ["R388","RP4","pCU1","R6K"]
# Step 1: Combine all dictionaries into a DataFrame
heatmap_data = np.zeros((16,24))
for i, column_to_sample_dicts in enumerate(all_dicts):
    for dict_idx, column_to_sample in enumerate(column_to_sample_dicts):
        for index_pair, day_label in column_to_sample.items():
            first, second = map(int, index_pair.split("-"))  # Split the index pair into two numbers
            heatmap_data[first-1,second-1,]=(i//2)+1
            # heatmap_data.append({
            #     "First": first,
            #     "Second": second,
            #     "Day": day_label,
            #     "Plasmid": "%s" % Plasmids[i // 2],
            #     "Dictionary": f"Dict {dict_idx + 1}"
            # })

# Step 5: Plot the heatmap
fig,ax = plt.subplots(1,1,figsize=(10, 6))
sns.heatmap(
    heatmap_data,
    annot=True,
    yticklabels=np.arange(1,17,1),
    xticklabels=np.arange(1,25,1),
    cmap="tab10",
    cbar=False,
    linewidths=0.5,
    ax=ax
)

# Add title and labels
ax.set_title("500-fold dilution", fontsize=16)
ax.set_xlabel("i7")
ax.set_ylabel("i5")

# Show the plot
fig.tight_layout()

"""
1:100 dilution
"""
# Keio #1 + R388
# three replicates; two conditions: with or without antibiotic pulse
R388_100_no_pulse = [
    {
        "1-1": "Day 0",
        "4-3": "Day 2",
        "1-10": "Day 3",
        "1-16": "Day 7",
        "1-22": "Day 10",
        "9-4": "Day 15",
        "9-10": "Day 20",
    },
    {
        "1-1": "Day 0",
        "5-3": "Day 2",
        "1-11": "Day 3",
        "1-17": "Day 7",
        "1-23": "Day 10",
        "9-5": "Day 15",
        "9-11": "Day 20",
    },
    {
        "1-1": "Day 0",
        "6-3": "Day 2",
        "1-12": "Day 3",
        "1-18": "Day 7",
        "1-24": "Day 10",
        "9-6": "Day 15",
        "9-12": "Day 20",
    },
]
R388_100_pulse = [
    {
        "1-1": "Day 0",
        "4-3": "Day 2",
        "5-10": "Day 3",
        "5-16": "Day 7",
        "5-22": "Day 10",
        "13-4": "Day 15",
        "13-10": "Day 20",
    },
    {
        "1-1": "Day 0",
        "5-3": "Day 2",
        "5-11": "Day 3",
        "5-17": "Day 7",
        "5-23": "Day 10",
        "13-5": "Day 15",
        "13-11": "Day 20",
    },
    {
        "1-1": "Day 0",
        "6-3": "Day 2",
        "5-12": "Day 3",
        "5-18": "Day 7",
        "5-24": "Day 10",
        "13-6": "Day 15",
        "13-12": "Day 20",
    },
]
R388_dicts_100 = [R388_100_no_pulse,R388_100_pulse]
save_data(data,R388_dicts_100,"R388_100_dilution",save=save)
# Keio #2 + RP4
RP4_100_no_pulse = [
    {
        "2-1": "Day 0",
        "4-4": "Day 2",
        "2-10": "Day 3",
        "2-16": "Day 7",
        "2-22": "Day 10",
        "10-4": "Day 15",
        "10-10": "Day 20",
    },
    {
        "2-1": "Day 0",
        "5-4": "Day 2",
        "2-11": "Day 3",
        "2-17": "Day 7",
        "2-23": "Day 10",
        "10-5": "Day 15",
        "10-11": "Day 20",
    },
    {
        "2-1": "Day 0",
        "6-4": "Day 2",
        "2-12": "Day 3",
        "2-18": "Day 7",
        "2-24": "Day 10",
        "10-6": "Day 15",
        "10-12": "Day 20",
    },
]
RP4_100_pulse = [
    {
        "2-1": "Day 0",
        "4-4": "Day 2",
        "6-10": "Day 3",
        "6-16": "Day 7",
        "6-22": "Day 10",
        "14-4": "Day 15",
        "14-10": "Day 20",
    },
    {
        "2-1": "Day 0",
        "5-4": "Day 2",
        "6-11": "Day 3",
        "6-17": "Day 7",
        "6-23": "Day 10",
        "14-5": "Day 15",
        "14-11": "Day 20",
    },
    {
        "2-1": "Day 0",
        "6-4": "Day 2",
        "6-12": "Day 3",
        "6-18": "Day 7",
        "6-24": "Day 10",
        "14-6": "Day 15",
        "14-12": "Day 20",
    },
]
RP4_dicts_100 = [RP4_100_no_pulse,RP4_100_pulse]
save_data(data,RP4_dicts_100,"RP4_100_dilution",save=save)
# Keio #3 + pCU1
pCU1_100_no_pulse = [
    {
        "3-1": "Day 0",
        "4-5": "Day 2",
        "3-10": "Day 3",
        "3-16": "Day 7",
        "3-22": "Day 10",
        "11-4": "Day 15",
        "11-10": "Day 20",
    },
    {
        "3-1": "Day 0",
        "5-5": "Day 2",
        "3-11": "Day 3",
        "3-17": "Day 7",
        "3-23": "Day 10",
        "11-5": "Day 15",
        "11-11": "Day 20",
    },
    {
        "3-1": "Day 0",
        "6-5": "Day 2",
        "3-12": "Day 3",
        "3-18": "Day 7",
        "3-24": "Day 10",
        "11-6": "Day 15",
        "11-12": "Day 20",
    },
]
pCU1_100_pulse = [
    {
        "3-1": "Day 0",
        "4-5": "Day 2",
        "7-10": "Day 3",
        "7-16": "Day 7",
        "7-22": "Day 10",
        "15-4": "Day 15",
        "15-10": "Day 20",
    },
    {
        "3-1": "Day 0",
        "5-5": "Day 2",
        "7-11": "Day 3",
        "7-17": "Day 7",
        "7-23": "Day 10",
        "15-5": "Day 15",
        "15-11": "Day 20",
    },
    {
        "3-1": "Day 0",
        "6-5": "Day 2",
        "7-12": "Day 3",
        "7-18": "Day 7",
        "7-24": "Day 10",
        "15-6": "Day 15",
        "15-12": "Day 20",
    },
]
pCU1_dicts_100 = [pCU1_100_no_pulse,pCU1_100_pulse]
save_data(data,pCU1_dicts_100,"pCU1_100_dilution",save=save)
# Keio #4 + R6K
R6K_100_no_pulse = [
    {
        "4-1": "Day 0",
        "4-6": "Day 2",
        "4-10": "Day 3",
        "4-16": "Day 7",
        "4-22": "Day 10",
        "12-4": "Day 15",
        "12-10": "Day 20",
    },
    {
        "4-1": "Day 0",
        "5-6": "Day 2",
        "4-11": "Day 3",
        "4-17": "Day 7",
        "4-23": "Day 10",
        "12-5": "Day 15",
        "12-11": "Day 20",
    },
    {
        "4-1": "Day 0",
        "6-6": "Day 2",
        "4-12": "Day 3",
        "4-18": "Day 7",
        "4-24": "Day 10",
        "12-6": "Day 15",
        "12-12": "Day 20",
    },
]
R6K_100_pulse = [
    {
        "4-1": "Day 0",
        "4-6": "Day 2",
        "8-10": "Day 3",
        "8-16": "Day 7",
        "8-22": "Day 10",
        "16-4": "Day 15",
        "16-10": "Day 20",
    },
    {
        "4-1": "Day 0",
        "5-6": "Day 2",
        "8-11": "Day 3",
        "8-17": "Day 7",
        "8-23": "Day 10",
        "16-5": "Day 15",
        "16-11": "Day 20",
    },
    {
        "4-1": "Day 0",
        "6-6": "Day 2",
        "8-12": "Day 3",
        "8-18": "Day 7",
        "8-24": "Day 10",
        "16-6": "Day 15",
        "16-12": "Day 20",
    },
]
R6K_dicts_100 = [R6K_100_no_pulse,R6K_100_pulse]
save_data(data,R6K_dicts_100,"R6K_100_dilution",save=save)

# confirm the dictionaries are correct
all_dicts = R388_dicts_100+RP4_dicts_100+pCU1_dicts_100+R6K_dicts_100
Plasmids = ["R388","RP4","pCU1","R6K"]
# Combine all dictionaries into a DataFrame
heatmap_data_100 = np.zeros((16,24))
for i, column_to_sample_dicts in enumerate(all_dicts):
    for dict_idx, column_to_sample in enumerate(column_to_sample_dicts):
        for index_pair, day_label in column_to_sample.items():
            first, second = map(int, index_pair.split("-"))  # Split the index pair into two numbers
            heatmap_data_100[first-1,second-1,]=(i//2)+1
# Plot the heatmap
fig2,ax2 = plt.subplots(1,1,figsize=(10, 6))
sns.heatmap(
    heatmap_data_100,
    annot=True,
    yticklabels=np.arange(1,17,1),
    xticklabels=np.arange(1,25,1),
    cmap="tab10",
    cbar=False,
    linewidths=0.5,
    ax=ax2
)

# Add title and labels
ax2.set_title("100-fold dilution", fontsize=16)
ax2.set_xlabel("i7")
ax2.set_ylabel("i5")

# Show the plot
fig2.tight_layout()
plt.show()