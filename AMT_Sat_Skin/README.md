# Satellite skin temperature append to input file
This code appends satellite based skin temperatures to the input files in input_data. This firstly appends data provided by Werenfrid Wimmer for Sentinel 3A, and then appends the CCI-SST data.

# Instructions
1. Run ISAR_SAT_MATCH.py to produce a csv file of ISAR matches to the AMT in situ observations on AMT28 and 29. This produces both 1 hr and 2 hr variants and the satellite validation statistics.
2. Run ISAR_SAT_APPEND.py to append the 2 hr variant to the initial input text files. This adds the '_S3A_ADDED' to the filename.
3. Run CCISST_run_ISAR_netcdf.py to append CCI-SST data to the ISAR netcdf provided by Werenfrid Wimmer. (This is needed for some plotting, but not the vertical temperature gradients analysis)
4. Run CCISST_run_AMT.py to append the CCI-SST to the input_data with the "_S3A_ADDED" component.

These steps get the input file from the raw data provided by Ming-xi Yang and Tom Bell, to the satellite skin temperature appended versions which are used in analysis...

# Optional Instructions for skin temperature validation...
