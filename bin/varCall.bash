#!/bin/bash
#================================================================
# varCall
#================================================================
#% SYNOPSIS
#+    varCall.bash ref read_1 read_2 mapped
#%
#% DESCRIPTION
#%    Determines where the sample differs from the reference genome
#%
#% INPUTS
#%    ref             path to input reference fasta file
#%    bam             path to first input fastq read file
#%    vcf             path to output vcf.gz file

# Error handling
set -eo pipefail

# Determines where the sample differs from the reference genome

ref=$1
bam=$2
vcf=$3
MAP_QUAL=$4
BASE_QUAL=$5
PLOIDY=$6

samtools index $bam
bcftools mpileup -q $MAP_QUAL -Q $BASE_QUAL -a INFO/AD -Ou -f $ref "$bam" |
    bcftools call --ploidy $PLOIDY -mf GQ - -Ou |
    bcftools norm -f $ref - -Oz -o "$vcf"
bcftools index $vcf
