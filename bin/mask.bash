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
#%    masked       output path to masked bed file INCLUDING rpt mask
#%    filtered     output path to masked bed file EXCLUDING rpt mask
#%    regions      output path to regions file (sites to keep)
#%    allsites     input path to allsites bed file

# Error handling
set -eo pipefail

# Inputs
rpt_mask=$1
vcf=$2
masked=$3
filtered=$4
regions_masked=$5
regions_filtered=$6
bam=$7
allsites=$8
MIN_READ_DEPTH=$9
MIN_ALLELE_FREQUENCY_ALT=${10}
MIN_ALLELE_FREQUENCY_REF=${11}

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

# exclude rpt mask
cat quality-mask.bed zerocov.bed | 
sort -k1,1 -k2,2n |
bedtools merge > $filtered

# Make bedfile of sites to keep - exc. mask positions
bedtools subtract -a $allsites -b $masked > $regions_masked

# Make bedfile of sites to keep - inc. mask positions
bedtools subtract -a $allsites -b $filtered > $regions_filtered

# Cleanup
rm quality-mask.vcf quality-mask.bed zerocov.bed
