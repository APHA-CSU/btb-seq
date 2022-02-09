#!/bin/bash

# Assigns cluster by matching patterns of cluster specific SNPs
# Compares SNPs identified in vcf file to lists in reference table

pair_id=$1
discrimPos=$2
stage1pat=$3
min_mean_cov=$4
min_cov_snp=$5
alt_prop_snp=$6
min_qual_snp=$7
min_qual_nonsnp=$8
pypath=$9
ref=$10

bcftools view -R ${discrimPos} -O v -o ${pair_id}.discrimPos.vcf ${pair_id}.vcf.gz
python3 $pypath/Stage1-test.py ${pair_id}_stats.csv ${stage1pat} $ref test ${min_mean_cov} ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} ${min_qual_nonsnp} ${pair_id}.discrimPos.vcf
mv _stage1.csv ${pair_id}_stage1.csv
