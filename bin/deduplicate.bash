#!/bin/bash
#================================================================
# Deduplicate
#================================================================
#% SYNOPSIS
#+    bov-tb reads_dir results_dir
#%
#% DESCRIPTION
#%    Removes potential duplicate data (sequencing and optical replcaites from the raw data set) from a pair of .fastq.gz files using FastUniq
#%
#% INPUTS
#%    pair_1          input path to first fastq.gz read file
#%    pair_2          input path to second fastq.gz read file
#%    deduplicated_1  output path to first deduplicated file
#%    deduplicated_2  output path to second deduplicated file

# Inputs
pair_1=$1
pair_2=$2
deduplicated_1=$3
deduplicated_2=$4

nl=$'\n'

# Unzip
gunzip -c $pair_1 > R1.fastq 
gunzip -c $pair_2 > R2.fastq

# Input List
echo "R1.fastq${nl}R2.fastq" > fqin.lst

# FastUniq
fastuniq -i fqin.lst -o $deduplicated_1 -p $deduplicated_2

# Cleanup
rm R1.fastq R2.fastq fqin.lst