#!/bin/bash
# Consensus calling and output snp table for snippy compatability 

pair_id=$1
ref=$2

 bcftools filter --SnpGap 5 -e 'DP<5 || INFO/AD[1]/(INFO/AD[1]+INFO/AD[0]) <0.8' ${pair_id}.norm.vcf.gz -Ob -o ${pair_id}.norm-flt.bcf
 bcftools index ${pair_id}.norm-flt.bcf
 bcftools consensus -f ${ref} -e 'TYPE="INDEL"' -m ${pair_id}_RptZeroMask.bed ${pair_id}.norm-flt.bcf |
  sed "/^>/ s/.*/>${pair_id}/" > ${pair_id}_consensus.fas
 echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > ${pair_id}_snps.tab
 bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' ${pair_id}.norm-flt.bcf |
 awk -F, '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$5":"$8+$9" "$4":"$6+$7}' >> ${pair_id}_snps.tab
