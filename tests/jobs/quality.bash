#!/bin/bash
#
#================================================================
# quality.bash
#================================================================
#
#% DESCRIPTION
#%    Checks the pipeline response against good and bad quality reads


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
python tests/utils/set_uniform_fastq_quality.py $quality /reads/$read1
python tests/utils/set_uniform_fastq_quality.py $quality /reads/$read2

# Zip 
gzip -f /reads/*.fastq

# Run nextflow
nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report "/artifacts/report.html"

# Check results
WGS_CLUSTER_CSV=$(sh tests/utils/print_todays_wgs_cluster.sh $name)
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV Outcome "$outcome"
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV flag "$flag"
python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV group "$group"
