#!/bin/bash

set -eo pipefail

# Args
ref=$1
read_1=$2
read_2=$3
mapped=$4

# Map to reference
bowtie2 -x $ref -1 $read_1 -2 $read_2 -p4 --end-to-end --no-mixed |
    samtools view -@2 -ShuF 2308 - |
    samtools sort -@2 - -o $mapped
