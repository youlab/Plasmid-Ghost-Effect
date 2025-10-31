import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
    with pd.ExcelWriter("%s.xlsx"%sample_name) as writer:
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
"""
1:100 dilution
"""
# three replicates; two conditions: with or without antibiotic pulse
R388_no_pulse = [
    {
        "9-13": "Day 0",
        "9-16": "Day 2",
        "10-13": "Day 3",
        "10-19": "Day 7",
        "11-13": "Day 10",
        "11-19": "Day 15",
        "12-13": "Day 20",
        "12-19": "Day 25",
        "13-13": "Day 30",
        "13-19": "Day 35",
    },
    {
        "9-13": "Day 0",
        "9-17": "Day 2",
        "10-14": "Day 3",
        "10-20": "Day 7",
        "11-14": "Day 10",
        "11-20": "Day 15",
        "12-14": "Day 20",
        "12-20": "Day 25",
        "13-14": "Day 30",
        "13-20": "Day 35",
    },
    {
        "9-13": "Day 0",
        "9-18": "Day 2",
        "10-15": "Day 3",
        "10-21": "Day 7",
        "11-15": "Day 10",
        "11-21": "Day 15",
        "12-15": "Day 20",
        "12-21": "Day 25",
        "13-15": "Day 30",
        "13-21": "Day 35",
    },
]
R388_pulse = [
    {
        "9-13": "Day 0",
        "9-16": "Day 2",
        "10-16": "Day 3",
        "10-22": "Day 7",
        "11-16": "Day 10",
        "11-22": "Day 15",
        "12-16": "Day 20",
        "12-22": "Day 25",
        "13-16": "Day 30",
        "13-22": "Day 35",
    },
    {
        "9-13": "Day 0",
        "9-17": "Day 2",
        "10-17": "Day 3",
        "10-23": "Day 7",
        "11-17": "Day 10",
        "11-23": "Day 15",
        "12-17": "Day 20",
        "12-23": "Day 25",
        "13-17": "Day 30",
        "13-23": "Day 35",
    },
    {
        "9-13": "Day 0",
        "9-18": "Day 2",
        "10-18": "Day 3",
        "10-24": "Day 7",
        "11-18": "Day 10",
        "11-24": "Day 15",
        "12-18": "Day 20",
        "12-24": "Day 25",
        "13-18": "Day 30",
        "13-24": "Day 35",
    },
]
R388_dicts = [R388_no_pulse,R388_pulse]
save_data(data,R388_dicts,"R388_100_dilution",save=save)

R388_dicts = R388_no_pulse+R388_pulse
# confirm the dictionaries are correct
heatmap_data = np.zeros((16,24))
for dict_idx, column_to_sample in enumerate(R388_dicts):
    for index_pair, day_label in column_to_sample.items():
        first, second = map(int, index_pair.split("-"))  # Split the index pair into two numbers
        heatmap_data[first-1,second-1,]=1
# Plot the heatmap
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
ax.set_title("100-fold dilution", fontsize=16)
ax.set_xlabel("i7")
ax.set_ylabel("i5")

# Show the plot
fig.tight_layout()
plt.show()