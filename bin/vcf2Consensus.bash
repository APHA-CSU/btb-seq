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

=======
# Inputs
ref=$1
bed=$2
vcf=$3
consensus=$4
snps=$5
MIN_READ_DEPTH=$6
MIN_ALLELE_FREQUENCY=$7
SNP_GAP=$8

# Filter
bcf=filtered.bcf
bcftools filter --SnpGap $SNP_GAP |
    -e "DP<${MIN_READ_DEPTH} || INFO/AD[1]/(INFO/AD[1]+INFO/AD[0]) < ${MIN_ALLELE_FREQUENCY} |
    || INFO/ADF == 0 || INFO/ADR == 0" |
    $vcf -Ob -o $bcf
bcftools index $bcf

# Call Consensus
base_name=`basename $consensus`
name="${base_name%%.*}"

bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $bed $bcf |
sed "/^>/ s/.*/>${name}/" > $consensus

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps

# Cleanup
rm $bcf $bcf.csi
