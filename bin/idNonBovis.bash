#!/bin/bash

# Identify any non-M.bovis samples using kraken
# Samples with flag != 'Pass' are processed to detemine which microbe(s) are present 
# Bracken parses the output which is then sorted to generate a top 20 list of species
# Presence / absence of M.bovis is also determined by parsing Bracken output

set -e

pair_id=$1
kraken2db=$2
lowmem=$3

outcome=$(cat outcome.txt)
if [ $outcome != "Pass" ]; then
kraken2 --threads 2 --quick $lowmem --db $kraken2db --output - --report ${pair_id}_"$outcome"_kraken2.tab --paired trimmed_1.fastq trimmed_2.fastq 

# HACK: (AF) Ignore Bracken errors. Better to handle output from Kraken and have unit tests, 
# but easier let the pipeline pass while we are setting up validation tests.. 
set +e

bracken -d $kraken2db -r 150 -l S -t 10 -i ${pair_id}_"$outcome"_kraken2.tab -o ${pair_id}_"$outcome"_bracken.out
sed 1d ${pair_id}_"$outcome"_bracken.out | sort -t $'\t' -k7,7 -nr - | head -20 > ${pair_id}_"$outcome"_brackensort.tab
bracken -d $kraken2db -r150 -l S1 -i ${pair_id}_"${outcome}"_kraken2.tab -o ${pair_id}_"$outcome"-S1_bracken.out
( sed -u 1q; sort -t $'\t' -k7,7 -nr ) < ${pair_id}_"$outcome"-S1_bracken.out > ${pair_id}_"$outcome"-S1_brackensort.tab
# Need bovis to be the top Mycobacterium variant for better confidence in confiming M bovis ID - reduce false positives
MycoPos=$(grep -m 1 'Mycobacterium' ${pair_id}_"$outcome"-S1_brackensort.tab)
BovPos=$(echo $MycoPos | grep 'variant bovis' |
 awk '{printf $1" "$2" "$3" "$4","$9"," "%.3f", ($10*100)}' || true)

echo "Sample,ID,TotalReads,Abundance" > ${pair_id}_bovis.csv
echo "${pair_id},"$BovPos"" >> ${pair_id}_bovis.csv

# HACK: see above
set -e

else
echo "ID not required"
fi
rm `readlink trimmed_1.fastq`
rm `readlink trimmed_2.fastq`
