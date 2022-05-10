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
#%    regions      input path to regions file containing all regions for masking/filtering(.bed)
#%    vcf          input path to variant (.vcf.gz)
#%    consensus    output path to consensus file (.fas)
#%    snps         output path to snps tab file (.tab)
#%    bcf          output path to filtered bcf file (.bcf)

# Error handling
set -eo pipefail

#=======
# Inputs
ref=$1
mask=$2
filter=$3
regions_masked=$4
regions_filtered=$5
vcf=$6
consensus=$7
snps=$8
bcf=$9
unmasked_consensus=${10}
MIN_ALLELE_FREQUENCY=${11}

# handle the case when the regions file is empty otherwise bcftools filter will faile
if [ ! -s $regions_masked ]; then
	# The file is empty.
	echo "LT708304-Mycobacteriumbovis-AF2122-97	-1	-1" > $regions_masked	
fi
if [ ! -s $regions_filtered ]; then
	# The file is empty.
	echo "LT708304-Mycobacteriumbovis-AF2122-97	-1	-1" > $regions_filtered	
fi

# Select SNPs
# applies filter/mask and chooses SNP sites
bcftools filter -R $regions_masked -i "ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) >= ${MIN_ALLELE_FREQUENCY}" $vcf -Ob -o $bcf
bcftools index $bcf

# applies filter and chooses SNP sites
bcftools filter -R $regions_filtered -i "ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) >= ${MIN_ALLELE_FREQUENCY}" $vcf -Ob -o unmasked_bcf
bcftools index unmasked_bcf

# Call Consensus
base_name=`basename $consensus`
name="${base_name%%.*}"

bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $mask $bcf |
sed "/^>/ s/.*/>${name}/" > $consensus

# TODO test unmasked_consensus
bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $filter unmasked_bcf |
sed "/^>/ s/.*/>${name}/" > $unmasked_consensus

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -R $regions_masked -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps
