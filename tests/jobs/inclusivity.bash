#!/bin/bash
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
accession=$(print_csv_value './tests/data/inclusivity_cases.csv' $TESTCASE accession)

# Fetch SRA Data
prefetch $accession -O ./
fasterq-dump ./$accession
rm ./$accession/*.sra
rm -r ./$accession

gzip ${accession}_1.fastq ${accession}_2.fastq
mv ${accession}_1.fastq.gz /reads/${accession}_S1_R1_001.fastq.gz
mv ${accession}_2.fastq.gz /reads/${accession}_S1_R2_001.fastq.gz

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "Pass"
assert_first_csv_row $WGS_CLUSTER_CSV "group" "$group"
