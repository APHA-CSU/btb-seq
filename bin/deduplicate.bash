#!/bin/bash

# This process removes potential duplicate data (sequencing and optical replcaites from the raw data set


pair_id=$1

nl=$'\n'

gunzip -c ${pair_id}_*_R1_*.fastq.gz > ${pair_id}_R1.fastq 
gunzip -c ${pair_id}_*_R2_*.fastq.gz > ${pair_id}_R2.fastq
echo "${pair_id}_R1.fastq${nl}${pair_id}_R2.fastq" > fqin.lst
fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
rm ${pair_id}_R1.fastq
rm ${pair_id}_R2.fastq