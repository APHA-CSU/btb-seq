
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

# Parameters
INDEL_GAP=5

# Inputs
ref=$1
bed=$2
vcf=$3
consensus=$4
snps=$5
VAR_QUAL=$6
MIN_READ_DEPTH=$7
MIN_ALT_PROPORTION=$8
window=$9

# Filter
bcf=filtered.bcf
masked_bcf=masked.bcf
filt_vcf=filtered.vcf
bcftools filter -e "TYPE!='snp' || %QUAL<${VAR_QUAL} || (AD[0]+AD[1])<=${MIN_READ_DEPTH} || ADF[1]==0 || ADR[1]==0 || (AD[1]/(AD[0]+AD[1]))<=${MIN_ALT_PROPORTION}" $vcf -Ob -o $bcf

vcftools --bcf $bcf --thin $window --recode-bcf --recode-INFO-all --out $masked_bcf

bcftools index $masked_bcf
bcftools view $masked_bcf -Oz -o $filt_vcf

# Call Consensus
base_name=`basename $consensus`
name="${base_name%%.*}"

bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $bed $masked_bcf |
sed "/^>/ s/.*/>${name}/" > $consensus

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $masked_bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps

# Cleanup
rm $bcf $masked_bcf $masked_bcf.csi

