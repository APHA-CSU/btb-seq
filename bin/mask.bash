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
vcf=$2
masked=$3
regions=$4
MIN_READ_DEPTH=$5
allsites=$6
MIN_ALLELE_FREQUENCY=$7

# Mask regions which don't have any evidence for ref or sufficient evidence for alt
bcftools filter -i "(ALT!='.' && INFO/AD[1] < ${MIN_READ_DEPTH}) ||
    (ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) < ${MIN_ALLELE_FREQUENCY}) ||
    (ALT='.' && AD=0)" $vcf -ov -o excluded-sites.vcf
bcftools filter -e 'ALT!="." && INFO/AD[0]/(INFO/AD[0]+INFO/AD[1]) > 0.6' excluded-sites.vcf -Ov -o quality-mask.vcf  
bedtools merge -i quality-mask.vcf > quality-mask.bed

## OLD FILTER # Find low coverage (<5 reads) regions
#bedtools genomecov -bga -ibam $bam |
#grep -w "[0-$((MIN_READ_DEPTH-1))]\$" |
#cat > low_cov.bed

# Merge with exisiting known repeat regions
cat quality-mask.bed $rpt_mask | 
sort -k1,1 -k2,2n |
bedtools merge > $masked

# Make bedfile of sites to keep
bedtools subtract -a $allsites -b $masked > $regions

# Cleanup
# rm excluded-sites.bed