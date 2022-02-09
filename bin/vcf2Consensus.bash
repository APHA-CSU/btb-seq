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
#%    mask         input path to mask (.bed)
#%    regions      input path to regions file containing all regions
#%                 for masking/filtering(.bed)
#%    vcf          input path to variant (.vcf.gz)
#%    consensus    output path to consensus file (.fas)
#%    snps         output path to snps tab file (.tab)
#%    bcf          output path to filtered bcf file (.bcf)

# Parameters
MIN_READ_DEPTH=5
MIN_ALLELE_FREQUENCY=0.8
INDEL_GAP=5

#=======
# Inputs
ref=$1
mask=$2
regions=$3
vcf=$4
consensus=$5
snps=$6
bcf=$7

# Filter
bcftools filter -R $regions $vcf -Ob -o $bcf
bcftools index $bcf

# Call Consensus
base_name=`basename $consensus`
name="${base_name%%.*}"

bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $mask $bcf |
sed "/^>/ s/.*/>${name}/" > $consensus

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -R $regions -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps
