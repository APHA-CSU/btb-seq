#!/bin/bash
#================================================================
# vcf2Consensus
#================================================================
#% SYNOPSIS
#+    mask.bash rpt_mask bam masked
#%
#% DESCRIPTION
#%    Consensus calling and output snp table for snippy compatability 
#%
#% INPUTS
#%    ref          input path to reference genome (.fas)
#%    bed          input path to mask (.bed)
#%    vcf          input path to variant (.vcf.gz)
#%    consensus    output path to consensus file (.fas)
#%    snps         output path to snps tab file (.tab)

set -e

# Parameters
MIN_READ_DEPTH=8
MIN_ALLELE_FREQUENCY=0.8
INDEL_GAP=5

#=======
# Inputs
ref=$1
bed=$2
vcf=$3
consensus=$4
snps=$5
bcf=$6

# Filter
bcftools filter --IndelGap $INDEL_GAP -i \
    "(DP>=${MIN_READ_DEPTH} && AD[1]>=${MIN_READ_DEPTH} && INFO/AD[1]/(INFO/AD[1]+INFO/AD[0]) >= ${MIN_ALLELE_FREQUENCY}) \
    || (DP>=${MIN_READ_DEPTH} && AD[0]>=${MIN_READ_DEPTH} && INFO/AD[0]/(INFO/AD[1]+INFO/AD[0]) >= ${MIN_ALLELE_FREQUENCY}) "\
    $vcf -Ob -o $bcf
bcftools index $bcf

# Bcf to vcf
bcftools view $bcf > filtered.vcf

# Add unfiltered output to mask
G=/home/aaronfishman/repos/btb-seq/references/full_genome.bed
F=filtered.vcf
M=$bed

bedtools subtract -a $G -b $F > subtracted.bed

# Add
cat subtracted.bed $M | 
sort -k1,1 -k2,2n |
bedtools merge > consensus_mask.bed

# Call Consensus
base_name=`basename $consensus`
name="${base_name%%.*}"

bcftools consensus -f ${ref} -e 'TYPE="indel"' -m consensus_mask.bed $bcf |
sed "/^>/ s/.*/>${name}/" > $consensus

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps
