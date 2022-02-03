#!/bin/bash #================================================================
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
ALT_MIN_READ_DEPTH=8
REF_MIN_READ_DEPTH=1
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
# AD[0] = (DP4[0]+DP4[1])
# AD[1] = (DP4[2]+DP4[3])


bcftools filter --IndelGap $INDEL_GAP -i \
    "(DP4[2]+DP4[3])>=${ALT_MIN_READ_DEPTH} && (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= ${MIN_ALLELE_FREQUENCY}" \
$vcf -Ob -o $bcf
bcftools index $bcf
bcftools view $bcf > snps.vcf

bcftools filter --IndelGap $INDEL_GAP -i \
    "(DP4[0]+DP4[1])>=${REF_MIN_READ_DEPTH} && (DP4[0]+DP4[1])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) >= ${MIN_ALLELE_FREQUENCY}" \
    $vcf -Ob -o refs.bcf
bcftools index refs.bcf
bcftools view refs.bcf > refs.vcf

# Concatenate
bcftools concat snps.vcf refs.vcf > filtered.vcf

# Full Genome Bes
printf "# Annotation file for storing repetitive regions with columns CHROM, START, END, RPT\n# number from 0, coordinates inclusive of start and stop\nLT708304-Mycobacteriumbovis-AF2122-97\t0\t4349905\t1" > full_genome.bed

#TODO: Better names
G=full_genome.bed
F=filtered.vcf
M=$bed

# Add unfiltered output to mask
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

bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' snps.vcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps
