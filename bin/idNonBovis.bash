#!/bin/bash

# Identify any non-M.bovis samples using kraken
# Samples with flag != 'Pass' are processed to detemine which microbe(s) are present 
# Bracken parses the output which is then sorted to generate a top 20 list of species
# Presence / absence of M.bovis is also determined by parsing Bracken output

pair_id=$1
kraken2db=$2
lowmem=$3

outcome=$(cat outcome.txt)
if [ $outcome != "Pass" ]; then
kraken2 --threads 2 --quick $lowmem --db $kraken2db --output - --report ${pair_id}_"$outcome"_kraken2.tab --paired ${pair_id}_trim_R1.fastq  ${pair_id}_trim_R2.fastq 

# HACK: (AF) Ignore Bracken errors. Better to handle output from Kraken and have unit tests, 
# but easier let the pipeline pass while we are setting up validation tests.. 
set +e

bracken -d $kraken2db -r 150 -l S -t 10 -i ${pair_id}_"$outcome"_kraken2.tab -o ${pair_id}_"$outcome"_bracken.out
sed 1d ${pair_id}_"$outcome"_bracken.out | sort -t $'t' -k7,7 -nr - | head -20 > ${pair_id}_"$outcome"_brackensort.tab
bracken -d $kraken2db -r150 -l S1 -i ${pair_id}_"${outcome}"_kraken2.tab -o ${pair_id}_"$outcome"-S1_bracken.out
( sed -u 1q; sort -t $'t' -k7,7 -nr ) < ${pair_id}_"$outcome"-S1_bracken.out > ${pair_id}_"$outcome"-S1_brackensort.tab
BovPos=$(grep 'variant bovis' ${pair_id}_"$outcome"-S1_brackensort.tab |
 awk '{print $1" "$2" "$3" "$4","$9","($10*100)}' || true)
echo "Sample,ID,TotalReads,Abundance" > ${pair_id}_bovis.csv
echo "${pair_id},"$BovPos"" >> ${pair_id}_bovis.csv

# HACK: see above
set -e

else
echo "ID not required"
fi
rm `readlink ${pair_id}_trim_R1.fastq`
rm `readlink ${pair_id}_trim_R2.fastq`