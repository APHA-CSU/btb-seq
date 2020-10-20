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

# Args
TESTCASE=$1

group=$(print_csv_value './tests/data/inclusivity_cases.csv' $TESTCASE group)
bovine_read1=$(print_csv_value './tests/data/inclusivity_cases.csv' $TESTCASE ftp_R1)
bovine_read2=$(print_csv_value './tests/data/inclusivity_cases.csv' $TESTCASE ftp_R2)


wget $bovine_read1 -O /reads/inclusivity-${group}_S1_R1_001.fastq.gz
wget $bovine_read2 -O /reads/inclusivity-${group}_S1_R2_001.fastq.gz

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "Pass"
assert_first_csv_row $WGS_CLUSTER_CSV "group" "$group"
