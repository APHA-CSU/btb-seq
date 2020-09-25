#!/bin/bash
#!/bin/sh
#
#================================================================
# tiny-reads-test.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset


# Set paths
mkdir -p /results/
mkdir -p /reads/tinyreads/
cp -r $PWD/tests/data/tinyreads/* /reads/tinyreads/

# Run nextflow
nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reads "/reads/MRP/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report "/artifacts/report.html"

# Check results
WGS_CLUSTER_CSV=/results/`sh tests/utils/print_todays_wgs_cluster.sh tinyreads`
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV Outcome InsufficientData
