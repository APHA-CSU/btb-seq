#!/bin/bash
#!/bin/sh
#
#================================================================
# quality.bash
#================================================================
#
#% DESCRIPTION
#%    Checks the pipeline response against good and bad quality reads

function test_quality {
    # Args
    quality=$1
    outcome=$2
    flag=$3
    group=$4

    # Set paths
    reads=/reads/$quality/
    mkdir -p $reads

    results=/results/$quality/
    mkdir -p $results

    # Download
    read1=B6-16_R1_26.fastq
    read2=B6-16_R2_26.fastq
    root=https://github.com/afishman/gitlfs/raw/master/lfs
    wget $root/$read1.gz
    wget $root/$read2.gz

    # Unzip
    gunzip -f ${read1}.gz ${read2}.gz

    # Set quality
    python tests/utils/set_uniform_fastq_quality.py $quality $read1 $read1
    python tests/utils/set_uniform_fastq_quality.py $quality $read2 $read2

    # Move over
    name=quality-${quality}
    mv $read1 $reads/${name}_S1_R1_001.fastq
    mv $read2 $reads/${name}_S1_R2_001.fastq

    # Zip 
    gzip -f $reads/*.fastq

    # Run nextflow
    nextflow run bTB-WGS_process.nf \
    --outdir "$results/" \
    --reads "$reads/*_{S*_R1,S*_R2}*.fastq.gz" \
    --lowmem '"--memory-map"' \
    -with-report "/artifacts/report-Q${quality}.html"
    
    # Check results
    WGS_CLUSTER_CSV=$results/`sh tests/utils/print_todays_wgs_cluster.sh $quality`
    python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV Outcome "$outcome"
    python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV flag "$flag"
    python tests/utils/assert_first_csv_row.py $WGS_CLUSTER_CSV group "$group"
}

# Check the two quality scores
test_quality 19 "CheckRequired" "LowCoverage" "NA"
test_quality 20 "Pass" "BritishbTB" "B6-16"
