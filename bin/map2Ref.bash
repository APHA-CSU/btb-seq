# map to reference sequence
# Aligns the individiual sequence reads to the reference genome 

pair_id=$1
ref=$2

bwa mem -M -t2 $ref  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq |
samtools view -@2 -ShuF 2308 - |
samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam