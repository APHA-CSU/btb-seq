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
read1=ERR4586795_1.fastq.gz
read2=ERR4586795_2.fastq.gz
root=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/005/ERR4586795/
wget $root/$read1 -O QT_S1_R1_001.fastq.gz -P /reads/
wget $root/$read2 -O QT_S1_R2_001.fastq.gz -P /reads/

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
