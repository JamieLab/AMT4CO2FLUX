# Vertical Temperature Gradients Analysis
This code produces the bulk air-sea CO2 fluxes using FluxEngine for AMT28 and AMT29, and then runs the analysis described in the supporting manuscript

# Instructions
1. Within folder AMT28 and AMT29 carry out these steps:
  1. Run FluxEngine_netcdfcreate.py - this converts the input text file into a netcdf and adds other parameters including the cool skin, corrected fCO2sw values etc. (produces AMT28_table_20min_version2_i.nc)
  2. Run FluxEngine_Ensemble_Generate.py - this produces a netcdf file where latitude is the data, and longitude is 100 ensembles of the data perturbed by the noise definitions. (Bodge way to do an ensemble analysis on 1d data with FluxEngine) This produces AMT28_table_20min_version2
  3. Run FluxEngine_run_and_save_data.py - this runs FluxEngine for all the possible combinations of vertical temperature gradients (no temp gradients, cool skin only, warm layer and cool skin), and saves it in FLUXENGINE_OUT.nc.

Now all the fluxes are calculated the averaging of these fluxes from 20 minutes to 3 hours is controlled in AMT_MERGED_PLOTTING.py
2. Run AMT_MERGED_PLOTTING.py to produce results (This requires the weight_stat.py file in the main directory, as well as AMT_mapo.py)
3. Run netcdf_data_packer.py to produce data outputs in Analysis_Output

# Optional Instructions for skin temperature validation...
