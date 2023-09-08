#!/bin/bash

#  readStats.bash
#
#  Created by ellisrichardj on 10/10/2019.
#  

# Error handling
set -eo pipefail

#  Define inputs - sample name
pair_id=$1

# Count reads in each catagory; in fastq files each read consists of four lines

    raw_R1=$(( $(zcat ${pair_id}_*_R1_*.fastq.gz | wc -l) / 4 )) # counts number of reads in file
    rm ${pair_id}_*_R1_*.fastq.gz
    rm ${pair_id}_*_R2_*.fastq.gz
    uniq_R1=$(( $(cat ${pair_id}_uniq_R1.fastq | wc -l) / 4 )) # counts number of reads in file
    rm `readlink ${pair_id}_uniq_R2.fastq`
    trim_R1=$(( $(cat ${pair_id}_trim_R1.fastq | wc -l) / 4 )) # counts number of reads in file
    decontam_R1=$(( $(zcat clean_R1.fastq.gz | wc -l) / 4 )) # counts number of reads in file
    num_map=$(samtools view -c ${pair_id}.mapped.sorted.bam) # samtools counts the number of mapped reads
    samtools depth -a ${pair_id}.mapped.sorted.bam > depth.txt # samtools outputs depth at each position
    avg_depth=$(awk '{sum+=$3} END { printf "%.3f", sum/NR}' depth.txt)
    zero_cov=$(awk 'BEGIN {count=0} $3<1 {++count} END {print count}' depth.txt)
    sites=$(awk '{++count} END {print count}' depth.txt)
    rm depth.txt
    


# Caluclate values and percentages

    num_raw=$(($raw_R1*2))
    num_uniq=$(($uniq_R1*2))
    num_trim=$(($trim_R1*2))
    num_clean=$(($decontam_R1*2))
    pc_aft_dedup=$(echo "scale=2; ($num_uniq*100/$num_raw)" |bc)
    pc_aft_trim=$(echo "scale=2; ($num_trim*100/$num_uniq)" |bc)
    pc_clean=$(echo "scale=2; ($num_clean*100/$num_trim)" |bc)
    pc_mapped=$(echo "scale=2; ($num_map*100/$num_trim)" |bc)
    genome_cov=$(echo "scale=2; (100-($zero_cov*100/$sites))" |bc)


# Define thresholds for flag assignment

    mindepth=10 # require an average of at least 10 reads per site 
    minpc=60 # require at least 60% of data maps to genome
    minreads=600000 # require at least 600,000 raw reads per sample
    minafttrim=60 # require at least 60% reads to pass quality filtering steps
    minclean=90 # require at least 90% of reads to be M.bovis

# This section assigns 'flags' based on the number of reads and the proportion mapping to reference genome
    
    if [[ ${pc_aft_trim%%.*} -lt $minafttrim ]]; then flag="LowQualData"
        elif [[ ${avg_depth%%.*} -ge $mindepth ]] && [[ ${pc_mapped%%.*} -gt $minpc ]]; then flag="Pass"
#        elif [[ ${pc_clean%%.*} -lt $minclean ]]; flag=""
        elif [[ ${avg_depth%%.*} -lt $mindepth ]] && [[ ${pc_mapped%%.*} -lt $minpc ]] && [[ $num_trim -gt $minreads ]]; then flag="Contaminated"
        elif [[ ${avg_depth%%.*} -lt $mindepth ]] && [[ $num_trim -lt $minreads ]]; then flag="InsufficientData"
#        elif [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="q_OtherMycobact"
        else flag="CheckRequired"
    fi


# Write values to csv file

    echo "Sample,NumRawReads,NumDedupReads,%afterDedup,NumTrimReads,%afterTrim,%Clean,NumMappedReads,%Mapped,MeanDepth,GenomeCov,Outcome" > ${pair_id}_stats.csv
    echo "${pair_id},"$num_raw","$num_uniq","$pc_aft_dedup","$num_trim","$pc_aft_trim","$pc_clean","$num_map","$pc_mapped","$avg_depth","$genome_cov","$flag"" >> ${pair_id}_stats.csv
    echo "$flag" > outcome.txt
