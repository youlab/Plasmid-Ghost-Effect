# Plasmid Ghost Effect
## Authors
Author: Zhengqing Zhou, Andrea Weiss, Zhixiang Yao, Kristen Lok, Hey-in Son, and Lingchong You

This is the github repo for reproducing analysis conducted in the paper *Dynamical memory underlies prolonged plasmid persistence after transient antibiotic treatment*

To ensure compatibility with the scripts in this repository, create the Conda environment using:

`conda env create -f py38_env.yml -n plasmid_ghost`

Raw 16S sequencing data are available under NCBI BioProject accession PRJNA1360409.

Below are the contents in each folder. 

## Folder: Simulations
Code for simulating and visualizing plasmid dynamics (Fig. 1 b, c and Fig. 4)

- ODE model details: `equations.py`

- clonal population: `clonal_theory.py`

- after antibiotic pulse: `pulse_effect.py`

- two-member comunity: `niche_partition.py`

## Folder: Clonal_experiments
Raw data (GFP/OD or selective plating) and code for processing, analyzing, and visualizing plasmid dynamics in clonal populations. (Fig. 1, 2, S1 - 3)

### Non-mobilizable plasmid dynamics
- `GFP_Plasmid_Processing.py`: process the raw data `./Raw_data/GE_LT5_GFP_calibration.xlsx` and `./Raw_data/GE_LT5_GFP.xlsx` through calibration of GFP/OD to plasmid abundance, and transform the dataset to `./LT_Data_py/*_mean.npy` and `./LT_Data_py/*_std.npy` files.
- `GFP_Plasmid_Dynamics.py`: subsequent time series visualization, quantification and visualization of plasmid half-lives based on the `.npy` files.

### Conjugative plasmid dynamics
- `Conjugative_Plating_Processing.py`: process raw data (selective plating) `./Raw_data/GE_LT10_plating.xlsx` that transforms into plasmid abundances with standard errors saved as `./LT_Data_py/*_mean.npy` and `./LT_Data_py/*_se.npy` files.
- `Conjugative_Plasmid_Dynamics.py`: time series visualization, quantification and visualization of plasmid half-lives.

### Non-mobilizable plasmid dynamics after antibiotic treatment
- `Does_Response_Calibration.py`: takes `./Raw_data/GE_LT18_GFP_calibration.xlsx` generate calibration curves between GFP/OD and plasmid abundances.
- `Dose_Response_processing.py`: takes `./Raw_data/GE_LT18_GFP.xlsx` to transform experimental GFP/OD readouts to plasmid abundances, with data saved as `./LT_Data_py/*_dose_response_mean.npy` and `./LT_Data_py/*_dose_response_std.npy` files.
- `Dose_Response_timeseries.py`: visualize selected time courses of plasmid abundance
- `Dose_Response_timeseries_all.py`: visualize all time courses
- `Dose_Response_visualization.py`: visualize the dose response of plasmid half-lives with antibiotic concentrations.

## Folder: Keio_Community
Includes raw data from selective plating of plasmid abundance; and barcode counts from next-generation amplicon sequencing and galaxy workflow processing (`Barcode_counting.ga`)
Includes scripts to process, analyze, and visualize the plasmid dynamics and community dynamics of synthetic *E. coli* Keio communities (Comm87 and Comm57). (Fig. 3, & S4 - 6)

### Comm87: 87-member community transferring plsamid R388
- `LT11_plating_process.py` processes the selective plating results from `./raw_data/GE_LT11_Keio.xlsx`, transforming them into `./processed_data/R388_mean.npy` and `./processed_data/R388_se.npy`.
- `LT11_plasmid_dynamics.py` takes `./processed_data/R388_mean.npy` and `./processed_data/R388_se.npy` to visualize plasmid dynamics and half-lives.
- `Merge_Count.py` and `organize_data.py` processes the barcode counts data under `raw_data` into `./processed_data/merged_counts.csv` and further process them as `./processed_data/R388_100_dilution.xlsx`. This excel file was further manually processed to exlude the strains not included in the experiment, and total reads and unmatched reads (from PhiX) as `./processed_data/R388_100_dilution_exclude_zeros.xlsx`.
- `Community_Dynamics_LT11.py` takes `./processed_data/R388_100_dilution_exclude_zeros.xlsx` to visualize strain-level community dynamics.

### Comm57: 57-member communities trasnferring one of the four plasmids: R388, RP4, pCU1, and R6K
- `LT21_plating_process.py` processes the selective plating results from `./raw_data/GE_LT21_Keio.xlsx`, transforming them into `./processed_data/*_mean.npy` and `./processed_data/*_se.npy`, of both 500-dilution and 100-dilution communities.
- `LT21_plasmid_dynamics_*.py` takes `./processed_data/*_mean.npy` and `./processed_data/*_se.npy` to visualize plasmid dynamics and half-lives.
- `Merge_Count.py` and `organize_data.py` processes the barcode counts data under `raw_data` into `./processed_data/merged_counts.csv` and further process them as `./processed_data/organized_NGS/*_dilution.xlsx`. 
- `Community_Dynamics_LT21_*.py` takes `./processed_data/organized_NGS/*_dilution.xlsx` to visualize strain-level community dynamics.

### Diveristy and donor abundances changes
The community dynamics datasets `Comm*.npy` were copied under the `coposition_all` folder, which were further used for community dynamics analyses.
- `diversity_change.py` takes the above datasets, and calculate the inverse Simpson index of the communities across time points.
- `donor_abundance_change.py` calculates the change in donor abundance between day 2 and 3.

## Folder: Sink_Community
16S sequencing and bioinformatic analyses were performed by SeqCenter LLC. This code base provides downstream visualization and analyses based on `./raw_data/16S_composition.xlsx` and GFP/OD readout in `./raw_data/Synk_LT1.xlsx`. Note that 16S compositions are only available for experiments with sponges, while GFP/OD data are available for experiments done with and without sponges. (Fig. 4 & S7)
- `16S_processing.py` takes in `./raw_data/16S_composition.xlsx`, cutoff to genus-evel composition, and organize the data to `./processed_data/*.npy`.
- `16S_composition.py` visualizes the community compositional changes.
- `LT1_processing.py` takes in GFP/OD readouts in `./raw_data/Synk_LT1.xlsx` and transform them into organized `*_mean.npy`and `*_std.npy` files.
- `LT1_plasmid_dynamics.py` visualizes plasmid dynamics over different conditions.

## Folder: Chemical_Treatment
Plasmid pSC101 dynamics in cloncal E. coli MG1655 under various chemical treatments. (Fig. S8)
- `Plating_Data_Processing.py` takes in selective plating results from `./Chemical_LT14.xlsx` and transform them into organized `*_mean.npy`and `*_se.npy` files.
- `Chemical_Curing_Plasmid_dynamics.py` visualizes plasmid dynamics and their half-lives over different conditions.

## License

This project is licensed under the MIT License.

## Contact

For questions, please contact:

- Lingchong You (lingchong.you@duke.edu)
- Department of Biomedical Engineering, Duke University
