# Determines where the sample differs from the reference genome

pair_id=$1
ref=$2

samtools index ${pair_id}.mapped.sorted.bam
bcftools mpileup -Q 10 -Ou -f $ref "${pair_id}.mapped.sorted.bam" |
bcftools call -mf GQ - -Ou |
bcftools norm -f $ref - -Oz -o "${pair_id}.norm.vcf.gz"