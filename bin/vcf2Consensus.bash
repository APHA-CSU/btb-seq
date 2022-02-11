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

# Inputs
ref=$1
mask=$2
regions=$3
vcf=$4
consensus=$5
snps=$6
MIN_READ_DEPTH=$7
MIN_ALLELE_FREQUENCY=$8
SNP_GAP=$9
SNP_WINDOW=${10}

# Filter
bcf=filtered.bcf
bcftools filter -R $regions --SnpGap $SNP_GAP -e 'ALT!="." && INFO/AD[0]/(INFO/AD[0]+INFO/AD[1]) > 0.6' $vcf |
#    -e "DP<${MIN_READ_DEPTH} ||
#    INFO/AD[1]/(INFO/AD[1]+INFO/AD[0]) < ${MIN_ALLELE_FREQUENCY} ||
#    INFO/ADF[1] == 0 || INFO/ADR[1] == 0" $vcf |
    bcftools +prune -w ${SNP_WINDOW}bp -n 1 -Ob -o $bcf
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

# Cleanup
rm $bcf $bcf.csi
