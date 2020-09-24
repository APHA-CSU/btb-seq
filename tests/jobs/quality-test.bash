#!/bin/bash
#!/bin/sh
#
#================================================================
# quality-test.bash
#================================================================
#
#% DESCRIPTION
#%    Checks the pipeline response against bad and good quality reads

function test_quality {
    quality=$1
    field=$2
    value=$3
    
    # Set paths
    reads=/reads/$quality/
    mkdir -p $reads

    results=/results/$quality/
    mkdir -p /results/$quality/

    # Download
    read1=B6-16_R1_26.fastq
    read2=B6-16_R2_26.fastq
    wget https://github.com/afishman/gitlfs/raw/master/lfs/$read1.gz
    wget https://github.com/afishman/gitlfs/raw/master/lfs/$read2.gz

    # Unzip
    gunzip -f ${read1}.gz ${read2}.gz

    # Change quality
    echo pwd $PWD
    ls tests/utils/
    python tests/utils/set_uniform_fastq_quality.py $quality $read1 $read1
    python tests/utils/set_uniform_fastq_quality.py $quality $read2 $read2

    # Move over
    mv $read1 $reads/quality-${quality}_S1_R1_001.fastq
    mv $read2 $reads/quality-${quality}_S1_R2_001.fastq

    # Zip 
    gzip -f $reads/*

    # Run nextflow
    nextflow run bTB-WGS_process.nf \
    --outdir "$results/" \
    --reads "$reads/*_{S*_R1,S*_R2}*.fastq.gz" \
    --lowmem '"--memory-map"' \
    -with-report "/results/report-$quality.html"
    
    # Check results
    WGS_CLUSTER_CSV=/results/`sh tests/utils/print_todays_wgs_cluster.sh`
    python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV $field $value
}

# Check the two quality scores
test_quality 19 Outcome Pass
test_quality 20 Outcome CheckRequired
