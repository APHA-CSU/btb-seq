#!/bin/bash

# This process removes potential duplicate data (sequencing and optical replcaites from the raw data set


dependpath=$1
pair_id=$2

echo Pair id $2

gunzip -c ${pair_id}_*_R1_*.fastq.gz > ${pair_id}_R1.fastq 
gunzip -c ${pair_id}_*_R2_*.fastq.gz > ${pair_id}_R2.fastq
echo "${pair_id}_R1.fastq\n${pair_id}_R2.fastq" > fqin.lst
${dependpath}/FastUniq/source/fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
rm ${pair_id}_R1.fastq
rm ${pair_id}_R2.fastq