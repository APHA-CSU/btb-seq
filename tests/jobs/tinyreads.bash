#!/bin/bash
#
#================================================================
# tinyreads.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a tiny dataset. 
#%    Asserts the Outcome is Insufficient Data

# Copy Data
cp -r $PWD/tests/data/tinyreads/* /reads/

# Run nextflow
nextflow run bTB-WGS_process.nf \
--outdir "/results" \
--reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report "/artifacts/report.html"

# Check results
WGS_CLUSTER_CSV=$(sh tests/utils/print_todays_wgs_cluster.sh tinyreads)
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV \
    Outcome InsufficientData
