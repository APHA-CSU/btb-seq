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

    cd /BovTB-nf/
    
    # Set paths
    reads=/reads/$quality/
    mkdir -p $reads

    results=/results/$quality/
    mkdir -p /results/$quality/

    # Download

    # Unzip
    gunzip /reads/$quality/...  /reads/$quality/...

    # Change quality
    python tests/set_uniform_fastq_quality.py $quality X Y
    python tests/set_uniform_fastq_quality.py $quality X Y

    # Zip 
    gzip X
    gzip Y

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
