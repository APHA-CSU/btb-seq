#!/bin/bash
#
#================================================================
# quality.bash
#================================================================
#
#% DESCRIPTION
#%    Checks the pipeline response against good or bad quality reads
#%    Input argument should be either 'low' or 'adequate'

# Import
source tests/utils/aliases.bash

TESTCASE=$1
echo if

# Test Case
if [ "$TESTCASE" == "low" ]; then
    echo Low quality test selected

    quality="19"
    outcome="CheckRequired"
    flag="LowCoverage"
    group="NA"

elif [ "$TESTCASE" == "adequate" ]; then
    echo Adequate quality test selected

    quality="20"
    outcome="Pass"
    flag="BritishbTB"
    group="B6-16"

else
    echo Unknown testcase: $TESTCASE
    exit 1

fi

# Download
name="B6-16"
read1=B6-16_SX_R1_26.fastq.gz
read2=B6-16_SX_R2_26.fastq.gz
root=https://github.com/afishman/gitlfs/raw/master/lfs
wget $root/$read1 $root/$read2 -P /reads/

# Unzip
gunzip -f /reads/*

# Set quality
python tests/utils/set_uniform_fastq_quality.py $quality /reads/*

# Zip 
gzip -f /reads/*

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "$outcome"
assert_first_csv_row $WGS_CLUSTER_CSV "flag" "$flag"
assert_first_csv_row $WGS_CLUSTER_CSV "group" "$group"
