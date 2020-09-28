#!/bin/bash
#
#================================================================
# quality.bash
#================================================================
#
#% DESCRIPTION
#%    Checks the pipeline response against good and bad quality reads

# Import
source tests/utils/aliases.bash

TESTCASE=$1

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
read1=B6-16_R1_26.fastq
read2=B6-16_R2_26.fastq
root=https://github.com/afishman/gitlfs/raw/master/lfs
wget $root/$read1.gz $root/$read2.gz -P /reads/

# Unzip
gunzip -f /reads/*

# Set quality
python tests/utils/set_uniform_fastq_quality.py $quality /reads/*

# Zip 
gzip -f /reads/*

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(sh tests/utils/print_todays_wgs_cluster.sh $name)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "$outcome"
assert_first_csv_row $WGS_CLUSTER_CSV "flag" "$flag"
assert_first_csv_row $WGS_CLUSTER_CSV "group" "$group"
