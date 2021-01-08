# map to reference sequence
# Aligns the individiual sequence reads to the reference genome 

dependpath=$1
pair_id=$2
ref=$3

$dependpath/bwa/bwa mem -M -t2 $ref  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq |
$dependpath/samtools-1.10/samtools view -@2 -ShuF 2308 - |
$dependpath/samtools-1.10/samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam