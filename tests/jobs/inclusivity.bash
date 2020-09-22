#!/bin/bash
#!/bin/sh
#
#================================================================
# inclusivity.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset


# Set paths
cd /BovTB-nf/
mkdir /results/
ln -s $PWD/00090/* /reads

# Run nextflow
nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reas "/reads/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report "/results/report.html"

# Check results
WGS_CLUSTER_CSV=/results/`sh tests/utils/print_todays_wgs_cluster.sh`
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV Outcome Pass
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV group B6-16
