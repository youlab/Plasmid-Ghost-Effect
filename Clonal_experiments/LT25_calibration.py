import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
plasmid_index = {"pSC101": 0, "colE1": 1, "pUC": 2}


def infer_plasmid_from_sheet(sheet_name: str) -> str | None:
    """Infer plasmid name from sheet label (e.g., 'Day3 pSC101')."""
    sheet_lower = sheet_name.lower()
    for plasmid in plasmid_index:
        if plasmid.lower() in sheet_lower:
            return plasmid
    return None

# Linear model
def linear_model(x, m, b):
    return m * x + b

def calibrate_sheet(xlsx_path, sheet_name, ax):
    """
    Load one sheet from Excel, subtract blank OD readings, calibrate each plasmid column,
    print report, and return fit results.
    
    Parameters:
    -----------
    xlsx_path : str
        Path to the Excel file
    sheet_name : str
        Name of the sheet to process
        ax : matplotlib.axes.Axes
            Axis to plot scatter points and fitted lines
    
    Returns:
    --------
    dict : {col: [slope, intercept, pcov]}
    """

    mix_ratio = np.array([100,50,0])
    
    # Load OD and GFP data
    od_df = pd.read_excel(xlsx_path,
                          sheet_name=sheet_name,
                          header=None, skiprows=3, nrows=4, usecols="B:J")
    gfp_df = pd.read_excel(xlsx_path,
                           sheet_name=sheet_name,
                           header=None, skiprows=9, nrows=4, usecols="B:J")
    
    od_np = od_df.to_numpy()
    gfp_np = gfp_df.to_numpy()
    
    # Blank readings (last row)
    od_blank = np.mean(od_np[-1, 0:3])
    gfp_blank = np.mean(gfp_np[-1, 0:3])
    
    # Calibration wells (rows 0-3)
    od_corr = od_np[0:3, :] - od_blank
    gfp_corr = gfp_np[0:3, :] - gfp_blank
    
    # Ratio GFP/OD
    y = gfp_corr / od_corr
    
    # Perform linear regression for each clone
    fit_result = {}
    

    x_fit = np.linspace(0, 105, 200)
    plasmid_name = infer_plasmid_from_sheet(sheet_name)
    plasmid_color = colors[plasmid_index[plasmid_name]] if plasmid_name is not None else "0.35"

    for i in range(9):
        y_i = y[:, i]
        # Linear fit
        popt, pcov = curve_fit(linear_model, mix_ratio, y_i)
        slope, intercept = popt
        slope_se, intercept_se = np.sqrt(np.diag(pcov))
        
        # Store results
        fit_result[i] = [slope, intercept, pcov]
        
        # Compute stats for report
        y_pred = linear_model(mix_ratio, slope, intercept)
        r2 = r2_score(y_i, y_pred)

        # Plot datapoints and fitted line for this clone on the provided axis.
        ax.scatter(mix_ratio, y_i, s=20, color=plasmid_color, alpha=0.9)
        ax.plot(x_fit, linear_model(x_fit, slope, intercept), lw=1, color=plasmid_color, alpha=0.8)
        
        # Print report
        print(f"[{sheet_name}] clone {i+1} calibration")
        report = (
            f"  Slope            : {slope:.3f} ± {slope_se:.3f} (SE)\n"
            f"  Intercept        : {intercept:.3f} ± {intercept_se:.3f} (SE)\n"
            f"  R²               : {r2:.3f}"
        )
        print(report)
        print()

    ax.set_title(sheet_name)
    ax.set_xlabel("P%")
    ax.set_ylabel("GFP/OD")
    ax.set_xlim(-5, 105)
    
    return fit_result

# Main processing
if __name__ == "__main__":
    xlsx_path = './raw_data/LT25_calibration.xlsx'
    mix_ratio = np.array([100, 50, 0])  # Updated mix ratios
    
    # Get sheet names from Excel file
    xls = pd.ExcelFile(xlsx_path)
    sheet_names = xls.sheet_names
    print(f"Found {len(sheet_names)} sheets: {sheet_names}\n")
    
    # Process first six sheets and plot each sheet on its own axis.
    selected_sheets = sheet_names[:6]
    n_cols = 3
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8,6))
    axes = axes.flat

    for idx, sheet_name in enumerate(selected_sheets):
        print(f"{'='*60}")
        print(f"Processing sheet: {sheet_name}")
        print(f"{'='*60}\n")

        calibration = calibrate_sheet(xlsx_path, sheet_name, axes[idx])
        output_pkl = f'./LT_Data_py/LT25_{sheet_name}_calibration.pkl'
        with open(output_pkl, 'wb') as f:
            pickle.dump(calibration, f)

    fig.tight_layout()
    # fig.savefig('./figures/LT25_calibration.png', dpi=300)
plt.show()

