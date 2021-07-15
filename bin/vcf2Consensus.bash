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
#%    vcf          input path to vcf file
#%    ref          input path to reference genome fas
#%    consensus    output path to consensus fas file
#%    snps         output path to snps tab file
#%    name         used as a label in the output consensus file

# Parameters
MIN_READ_DEPTH=5
MIN_ALLELE_FREQUENCY=0.8
INDEL_GAP=5

# Inputs
ref=$1
bed=$2
vcf=$3
consensus=$4
snps=$5
name=$6

# Filter
bcf=filtered.bcf
bcftools filter --IndelGap $INDEL_GAP -e "DP<${MIN_READ_DEPTH} && AF<${MIN_ALLELE_FREQUENCY}" $vcf -Ob -o $bcf
bcftools index $bcf

# Call Consensus
bed=consensus.bed
bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $bed $bcf |
sed "/^>/ s/.*/>${name}/" > $consensus

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps

# Cleanup
rm $bcf $bcf.csi