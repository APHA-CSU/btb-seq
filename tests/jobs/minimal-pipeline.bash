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
cd /BovTB-nf/
mkdir /results/
ln -s $PWD/tests/data/minimal-read-pair/ /reads

# Run nextflow
nextflow run bTB-WGS_process.nf --outdir "/results/" --reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz"

# Check results
WGS_CLUSTER_CSV=/results/`sh tests/utils/print_todays_wgs_cluster.sh`
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV Outcome InsufficientData
