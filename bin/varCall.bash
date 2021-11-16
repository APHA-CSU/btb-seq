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


# Determines where the sample differs from the reference genome

ref=$1
bam=$2
vcf=$3

samtools index $bam
bcftools mpileup -Q 10 -Ou -f $ref "$bam" |
bcftools call -mf GQ - -Ou |
bcftools norm -f $ref - -Oz -o "$vcf"
