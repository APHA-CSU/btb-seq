# Masking known repeats regions and sites with zero coverage
# Ensure that consensus only includes regions of the genome where there is high confidence */

dependpath=$1
pair_id=$2
rptmask=$3

$dependpath/bedtools2/bin/bedtools genomecov -bga -ibam ${pair_id}.mapped.sorted.bam |
grep -w "0\$" > ${pair_id}_zerocov.bed
cat ${pair_id}_zerocov.bed $rptmask | sort -k1,1 -k2,2n |
$dependpath/bedtools2/bin/bedtools merge > ${pair_id}_RptZeroMask.bed