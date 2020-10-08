#!/bin/bash
#
#================================================================
# lod.bash
#================================================================
#
#% DESCRIPTION
#%    Checks the pipeline response against specimens computationally mixed 
#%    with real-world M. Bovis and Avium samples

# Import
source tests/utils/aliases.bash

# Params
total_num_reads=6.5e5
outcomes=("Pass" "Pass" "CheckRequired" "CheckRequired")
props=("1.0" "0.55" "0.45" "0.0")

# Args
TESTCASE=$1
outcome=${outcomes[$TESTCASE]}
prop=${props[$TESTCASE]}

# Download files
mkdir -p fastq
bovine_read1=ref-genome-sample_SX_R1_001.fastq.gz
bovine_read2=ref-genome-sample_SX_R2_001.fastq.gz

avium_read1=SRR10541896_1.fastq.gz
avium_read2=SRR10541896_2.fastq.gz

bovine_read1=SRR10541896_1.fastq.gz
bovine_read2=SRR10541896_2.fastq.gz

avium_root=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/096/SRR10541896/
bovine_root=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/096/SRR10541896/

wget -P fastq \
  $bovine_root/$bovine_read1 \
  $bovine_root/$bovine_read2 \
  $avium_root/$avium_read1 \
  $avium_root/$avium_read2

# TODO: Combine
echo combining..
name=lod-bovine-${prop}
combine_fastq --seed 1 \
    $total_num_reads \
    $prop \
    fastq/$bovine_read1 \
    fastq/$avium_read1 \
    /reads/${name}_SX_R1_X.fastq.gz

combine_fastq --seed 1 \
    $total_num_reads \
    $prop \
    fastq/$bovine_read2 \
    fastq/$avium_read2 \
    /reads/${name}_SX_R2_X.fastq.gz

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "$outcome"
