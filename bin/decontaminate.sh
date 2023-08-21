#!/bin/bash

set -eo pipefail


read_1=$2
read_2=$3
$ref=$3

bbduk.sh in1=$read_1 in2=$read_2 \
    outm1=clean_R1.fastq.gz outm2=clean_R2.fastq.gz \
    ref=$ref \
    stats=stats.txt mcf=0.5
