# Masking known repeats regions and sites with zero coverage
# Ensure that consensus only includes regions of the genome where there is high confidence */

pair_id=$1
rptmask=$2

bedtools genomecov -bga -ibam ${pair_id}.mapped.sorted.bam |
grep -w "0\$" | cat > ${pair_id}_zerocov.bed
cat ${pair_id}_zerocov.bed $rptmask | sort -k1,1 -k2,2n |
bedtools merge > ${pair_id}_RptZeroMask.bed
