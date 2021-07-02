#!/bin/bash
#================================================================
# Trim
#================================================================
#% SYNOPSIS
#+    trim.bash read_1 read_2 trimmed_1 trimmed_2
#%
#% DESCRIPTION
#%    Trim adapters and low quality bases from fastq data
#%    Removes the adapters which are added during the lab processing and and any low quality data
#%
#% INPUTS
#%    adapters        path to adapter sequence
#%    read_1          input path to first fastq read file
#%    read_2          input path to second fastq read file
#%    trimmed_1       output path to first trimmed fastq file
#%    trimmed_2       output path to second trimmed fafstq file

adapter=$1
read_1=$2
read_2=$3
trimmed_1=$4
trimmed_2=$5

# Error if adapter does not exist, as trimmomatic won't!
if [[ ! -f $adapter ]]; then
    echo $adapter does not exist
    exit 1
fi

java -jar /usr/local/bin/trimmomatic.jar PE \
    -threads 2 \
    -phred33 \
    $read_1 \
    $read_2  \
    $trimmed_1 \
    fail1.fastq \
    $trimmed_2 \
    fail2.fastq \
    ILLUMINACLIP:$adapter:2:30:10 \
    SLIDINGWINDOW:10:20 \
    MINLEN:36

rm fail1.fastq fail2.fastq