#!/bin/bash
#================================================================
# Mask
#================================================================
#% SYNOPSIS
#+    mask.bash rpt_mask bam masked
#%
#% DESCRIPTION
#%    Masking known repeats regions and sites with zero coverage
#%    Ensure that consensus only includes regions of the genome where there is high confidence
#%
#% INPUTS
#%    rpt_mask     path to repeat mask
#%    bam          input path to mapped bam file
#%    masked       output path to masked bed file

# Inputs
rpt_mask=$1
bam=$2
masked=$3

# Find zero coverage regions
bedtools genomecov -bga -ibam $bam |
grep -w "0\$" | 
cat > zero_cov.bed

# Mask repeat regions
cat zero_cov.bed $rpt_mask | 
sort -k1,1 -k2,2n |
bedtools merge > $masked

# Cleanup
rm zero_cov.bed