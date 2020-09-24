#!/bin/bash
#!/bin/sh
#
#================================================================
# minimal-pipeline.sh
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset


# Set paths
mkdir -p /results/
mkdir -p /reads/MRP/
cp -r $PWD/tests/data/minimal-read-pair/* /reads/MRP/

# Run nextflow
nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reads "/reads/MRP/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report "/results/reports/report.html"

# Check results
WGS_CLUSTER_CSV=/results/`sh tests/utils/print_todays_wgs_cluster.sh MRP`
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV Outcome InsufficientData
