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
#%    vcf          input path to vcf file
#%    masked       output path to masked bed file
#%    regions      output path to regions file (sites to keep)
#%    allsites     input path to allsites bed file

# Error handling
set -eo pipefail

# Inputs
rpt_mask=$1
vcf=$2
masked=$3
regions=$4
bam=$5
allsites=$6
MIN_READ_DEPTH=$7
MIN_ALLELE_FREQUENCY_ALT=$8
MIN_ALLELE_FREQUENCY_REF=$9

# Construct a mask: 
# mask regions which don't have {sufficient evidence for alt AND sufficient evidence for the REF}
bcftools filter -i "(ALT!='.' && INFO/AD[1] < ${MIN_READ_DEPTH} && INFO/AD[0]/(INFO/AD[0]+INFO/AD[1]) <= ${MIN_ALLELE_FREQUENCY_REF}) ||
    (ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) < ${MIN_ALLELE_FREQUENCY_ALT} && INFO/AD[0]/(INFO/AD[0]+INFO/AD[1]) <= ${MIN_ALLELE_FREQUENCY_REF})" $vcf -ov -o quality-mask.vcf
bedtools merge -i quality-mask.vcf > quality-mask.bed

# mash regions where there is zero coverage
bedtools genomecov -bga -ibam $bam | grep -w "0\$" | cut -f -3 > zerocov.bed

# Merge with exisiting known repeat regions
cat quality-mask.bed zerocov.bed $rpt_mask | 
sort -k1,1 -k2,2n |
bedtools merge > $masked

# Make bedfile of sites to keep
bedtools subtract -a $allsites -b $masked > $regions

# Cleanup
rm quality-mask.vcf quality-mask.bed
