# Plasmid Ghost Effect
## Authors
Author: Zhengqing Zhou, Andrea Weiss, Zhixiang Yao, Xiaoli Chen, Jing-Mei Qian, Kristen Lok, Grayson S. Hamrick, Hye-in Son, and Lingchong You

This is the github repo for reproducing analysis conducted in the paper *Ecological memory of prolonged plasmid persistence after transient antibiotic exposure*

Raw amplicon sequencing data for E. coli Keio communities are available on Zenodo https://zenodo.org/records/18226472.

Raw 16S sequencing data are available under NCBI BioProject accession PRJNA1360409.

Below are the contents in each folder. 

## Requirements
To ensure compatibility with the scripts in this repository, create the Conda environment using:

`conda env create -f environment.yml -n plasmid_ghost`

The scripts have been tested on:<br>
Python: 3.14.2<br>
numpy: 2.4.1<br>
scipy: 1.17.0<br>
pandas: 2.3.3<br>
matplotlib: 3.10.8<br>
scikit-learn: 1.8.0

## Folder: Simulations
Code for simulating and visualizing plasmid dynamics (Fig. 1 b, c and Fig. 4)

- ODE model details: `equations.py`

- clonal population: `clonal_theory.py`

- after antibiotic pulse: `pulse_effect.py`

- two-member comunity: `niche_partition.py`

## Folder: Clonal_experiments
Raw data (GFP/OD or selective plating) and code for processing, analyzing, and visualizing plasmid dynamics in clonal populations, including ancestors and evolved clones. (Fig. 1, 2, S1, S2, S5, S7)

### Non-mobilizable plasmid dynamics
- Generates Fig. 1 & S1
- `GFP_Plasmid_Processing.py`: process the raw data `./Raw_data/GE_LT5_GFP_calibration.xlsx` and `./Raw_data/GE_LT5_GFP.xlsx` through calibration of GFP/OD to plasmid abundance, and transform the dataset to `./LT_Data_py/*_mean.npy` and `./LT_Data_py/*_std.npy` files.
- `GFP_Plasmid_Dynamics.py`: subsequent time series visualization, quantification and visualization of plasmid half-lives based on the `.npy` files.

### Conjugative plasmid dynamics
- Generates Fig. 1 & S2
- `Conjugative_Plating_Processing.py`: process raw data (selective plating) `./Raw_data/GE_LT10_plating.xlsx` that transforms into plasmid abundances with standard errors saved as `./LT_Data_py/*_mean.npy` and `./LT_Data_py/*_se.npy` files.
- `Conjugative_Plasmid_Dynamics.py`: time series visualization, quantification and visualization of plasmid half-lives.

### Non-mobilizable plasmid dynamics after antibiotic treatment
- Generates Fig. 2 & S5
- `Does_Response_Calibration.py`: takes `./Raw_data/GE_LT18_GFP_calibration.xlsx` generate calibration curves between GFP/OD and plasmid abundances.
- `Dose_Response_processing.py`: takes `./Raw_data/GE_LT18_GFP.xlsx` to transform experimental GFP/OD readouts to plasmid abundances, with data saved as `./LT_Data_py/*_dose_response_mean.npy` and `./LT_Data_py/*_dose_response_std.npy` files.
- `Dose_Response_timeseries.py`: visualize selected time courses of plasmid abundance
- `Dose_Response_timeseries_all.py`: visualize all time courses
- `Dose_Response_visualization.py`: visualize the dose response of plasmid half-lives with antibiotic concentrations.

### Evolved clones dynamics from antibiotic treatment experiment
- Generates Fig. S7
- `LT25_calibration.py`: read day 0 data with GFP/OD for P_0 = 100%, 50%, and 0%. Generate a calibration curve based on these datapoints `LT25_Day*_calibration.pkl` (* being 3 / 20; pSC101, colE1, and pUC).
- `LT25_processing.py`: process the raw data `./Raw_data/LT25_*.xlsx` files (* being pSC101, colE1, and pUC) by loading the calibration curve `.pkl` files to generate mean and std plasmid abundances `./LT_Data_py/*_mean.npy` and `./LT_Data_py/*_std.npy`.
- `LT25_Dynamics.py`: subsequent time series visualization, quantification and visualization of plasmid half-lives based on the `.npy` files.

## Folder: qPCR
Raw data and code for analyzing qPCR data to quantify plasmid copy number. Chromosome copy was quantified by primers targeting dxs gene; pSC101 and colE1 copy was quantified by primers targeting GFP; pUC copy was quantified by primers targeting its backbone. In particular, pSC101 was a concatemer, bearing two copies of GFP per molecule, thus its copy number was divided by 2. (Fig. S3)

- `qPCR_process_*.py` first generates calibration curves, which contains the amplification efficiency information of the two targets (chromosome and plasmid); then load the actual samples to infer the copy number.
- `visualize_PCN.py` visualizes plasmid copy number across 3 biological replicates.

## Folder: Growth_data
Raw data and code for analyzing growth curves and plasmid burden. (Fig. S4 & S6)

- `growth_curve_extraction_*.py` files process raw growth curve data through blanking, and visualize them in the original 96-well layout.

### Ancestor plasmid burden
- Generates Fig. S4
- `GFP_ancestor_growth_curve.py` and `GFP_ancestor_fitness.py` visualize MG1655 and plasmid-carrying (pSC101, colE1, and pUC) strains growth curve and estimate the fitness cost of plasmid on MG1655.
- `conjugative_ancestor_growth.py` and `conjugative_ancestor_fitness.py` visualize DA28102 and plasmid-carrying (pCU1 and R6K) strains growth curve and estimate the fitness cost of plasmid on MG1655.

### Evolved clone burden
- Generates Fig. S6
- `GFP_evolved_fitness.py` visualizes growth curves of evolved clones on day 3 and day 20 of the antibiotic treatment experiment, and quantifies plasmid fitness effects.
- `fitness_over_time.py` visualizes fitness effects of ancestral and evolved clones.

## Folder: Keio_Community
Includes raw data from selective plating of plasmid abundance; and barcode counts from next-generation amplicon sequencing and galaxy workflow processing (`Barcode_counting.ga`)
Includes scripts to process, analyze, and visualize the plasmid dynamics and community dynamics of synthetic *E. coli* Keio communities (Comm87 and Comm57). (Fig. 3, 4, S8 - S12)

### Folder Comm87: 87-member community transferring plsamid R388
- `LT11_plating_process.py` processes the selective plating results of non-treated and Ab-pulsed communities from `./raw_data/GE_LT11_Keio.xlsx`, transforming them into `./processed_data/R388_mean.npy` and `./processed_data/R388_se.npy`.
- `LT22_plating_process.py` processes invaded community by plasmid-free strain 76 from `./raw_data/LT22_invasion_plating.xlsx`, transforming them into `./processed_data/R388_inv_mean.npy` and `./processed_data/R388_inv_se.npy`.
- `post_Ab_plasmid_dynamics.py` takes `./processed_data/R388_mean.npy` and `./processed_data/R388_se.npy` to visualize R388 dynamics and half-lives before and after antibiotic pulse.
- `invasion_plasmid_dynamics.py` takes `./processed_data/R388_mean.npy`, `./processed_data/R388_se.npy`, `./processed_data/R388_inv_mean.npy` and `./processed_data/R388_inv_se.npy`, to visualize plasmid dynamics after antibiotic pulse with or without the plasmid-free strain 76 invasion; and visualizes plasmid half-lives before and after antibiotic pulse, and after invasion.
- `Merge_Count.py` and `organize_data.py` processes the barcode counts data under `raw_data` into `./processed_data/merged_counts.csv` and further process them as `./processed_data/R388_100_dilution.xlsx`. This excel file was further manually processed to exlude the strains not included in the experiment, and total reads and unmatched reads (from PhiX) as `./processed_data/R388_100_dilution_exclude_zeros.xlsx`.
- `Community_Dynamics_LT11.py` takes `./processed_data/R388_100_dilution_exclude_zeros.xlsx` to visualize strain-level community dynamics.

### Comm57: 57-member communities trasnferring one of the four plasmids: R388, RP4, pCU1, and R6K
- `LT21_plating_process.py` processes the selective plating results from `./raw_data/GE_LT21_Keio.xlsx`, transforming them into `./processed_data/*_mean.npy` and `./processed_data/*_se.npy`, of both 500-dilution and 100-dilution communities.
- `LT22_invasion_plating_process.py` processes invaded community by plasmid-free strain 3 (into Comm57+pCU1) and 4 (into Comm57+R6K) from `./raw_data/LT22_invasion_plating.xlsx`, transforming them into `./processed_data/*_inv_mean.npy` and `./processed_data/*_inv_se.npy`. * being pCU1 or R6K.
- `LT21_postAb_plasmid_dynamics_*.py` takes `./processed_data/*_mean.npy` and `./processed_data/*_se.npy` to visualize plasmid dynamics and half-lives before and after antibiotic pulse. * being the dilution rate (100 or 500).
- `invasion_plasmid_dynamics.py` takes `./processed_data/*_mean.npy`, `./processed_data/*_se.npy`, `./processed_data/*_inv_mean.npy` and `./processed_data/*_inv_se.npy`, to visualize plasmid dynamics after antibiotic pulse with or without plasmid-free strain invasion; and visualizes plasmid half-lives before and after antibiotic pulse, and after invasion.
- `Merge_Count.py` and `organize_data.py` processes the barcode counts data under `raw_data` into `./processed_data/merged_counts.csv` and further process them as `./processed_data/organized_NGS/*_dilution.xlsx`. 
- `Community_Dynamics_*.py` takes `./processed_data/organized_NGS/*_dilution.xlsx` to visualize strain-level community dynamics.

### Diveristy and donor abundances changes
The community dynamics datasets `Comm*.npy` were copied under the `coposition_all` folder, which were further used for community dynamics analyses.
- `diversity_change.py` takes the above datasets, and calculate the inverse Simpson index of the communities across time points.
- `donor_abundance_change.py` calculates the change in donor abundance between day 2 and 3.

## Folder: Sink_Community
16S sequencing and bioinformatic analyses were performed by SeqCenter LLC. This code base provides downstream visualization and analyses based on `./raw_data/16S_composition.xlsx` and GFP/OD readout in `./raw_data/Synk_LT1.xlsx`. Note that 16S compositions are only available for experiments with sponges, while GFP/OD data are available for experiments done with and without sponges. (Fig. 5 & S13)
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
