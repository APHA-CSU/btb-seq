#!/bin/bash
#
#================================================================
# aliases.bash
#================================================================
#
#% DESCRIPTION
#%    A number of aliases useful for testing
shopt -s expand_aliases

alias nextflowtest="bash tests/utils/nextflowtest.bash"

alias print_todays_wgs_cluster="sh tests/utils/print_todays_wgs_cluster.sh"

alias assert_first_csv_row="python tests/utils/assert_first_csv_row.py"

alias combine_fastq="python tests/utils/combine_fastq.py"

alias print_csv_value="python tests/utils/print_csv_value.py"