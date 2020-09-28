#!/bin/bash
#
#================================================================
# aliases.bash
#================================================================
#
#% DESCRIPTION
#%    A number of aliases useful for testing


alias nextflowtest="nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report /artifacts/report.html"

alias print_todays_wgs_cluster="sh tests/utils/print_todays_wgs_cluster.sh"

alias assert_first_csv_row="python tests/utils/assert_first_csv_row.py"
