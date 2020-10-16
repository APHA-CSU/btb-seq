#!/bin/bash
#!/bin/sh
#
#================================================================
# inclusivity.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset


# Import
source tests/utils/aliases.bash

bovine_root=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/005/ERR4586795/
bovine_read1=ERR4586795_1.fastq.gz
bovine_read2=ERR4586795_2.fastq.gz

wget $bovine_root/$bovine_read1 -O /reads/inclusivity-test_S1_R1_001.fastq.gz
wget $bovine_root/$bovine_read2 -O /reads/inclusivity-test_S1_R2_001.fastq.gz

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "Pass"
assert_first_csv_row $WGS_CLUSTER_CSV "group" "B6-16"

