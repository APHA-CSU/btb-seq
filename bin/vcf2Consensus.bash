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
regions=$3
vcf=$4
consensus=$5
snps=$6
bcf=$7
MIN_ALLELE_FREQUENCY=$8
pair_id=$9
publishDir=${10}


# handle the case when the regions file is empty otherwise bcftools filter will fail
if [ ! -s $regions ]; then
	# The file is empty.
	echo "LT708304-Mycobacteriumbovis-AF2122-97	-1	-1" > $regions	
fi

# Select SNPs
# applies filter/mask and chooses SNP sites
bcftools filter -i "ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) >= ${MIN_ALLELE_FREQUENCY}" $vcf -Ob -o $bcf
bcftools index $bcf

# Call Consensus
bcftools consensus -f ${ref} -e 'TYPE="indel"' -m $mask $bcf |
sed "/^>/ s/.*/>${pair_id}/" > $consensus

# Count Ns in consensus file and extract submission number from sample name
ncount=$(grep -o 'N' $consensus | wc -l)
set +e; submission_number=$(echo $pair_id | grep -Eo '[0-9]{2,2}\-[0-9]{4,5}\-[0-9]{2,2}'); set -e
echo -e "Sample,Submission,Ncount,ResultLoc" > ${pair_id}_ncount.csv
if [ -n "$submission_number" ]; then
    echo -e "${pair_id},$submission_number,$ncount,$publishDir" >> ${pair_id}_ncount.csv
else
    echo -e "${pair_id},${pair_id},$ncount,$publishDir" >> ${pair_id}_ncount.csv
fi

# Write SNPs table
echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > $snps

bcftools query -R $regions -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' $bcf |
awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> $snps
