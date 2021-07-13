#!/bin/bash
#================================================================
# map2Ref
#================================================================
#% SYNOPSIS
#+    map2Ref.bash ref read_1 read_2 mapped
#%
#% DESCRIPTION
#%    Map to reference sequence. Aligns the individiual sequence reads to the reference genome
#%
#% INPUTS
#%    ref             path to input reference fasta file
#%    read_1          path to first input fastq read file
#%    read_2          path to second input fastq read file
#%    mapped          path to output mapped bam file

# Args
ref=$1
read_1=$2
read_2=$3
mapped=$4

# Map to reference
bwa mem -M -t2 $ref $read_1 $read_2 |
samtools view -@2 -ShuF 2308 - |
samtools sort -@2 - -o $mapped