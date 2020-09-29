#!/bin/bash
#
#================================================================
# tinyreads.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a tiny dataset. 
#%    Asserts the Outcome is Insufficient Data

# Import
source tests/utils/aliases.bash

# Copy Data
cp -r $PWD/tests/data/tinyreads/* /reads/

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV \
    Outcome InsufficientData
