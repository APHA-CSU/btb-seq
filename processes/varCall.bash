# Determines where the sample differs from the reference genome

dependpath=$1
pair_id=$2

$dependpath/samtools-1.10/samtools index ${pair_id}.mapped.sorted.bam
bcftools mpileup -Q 10 -Ou -f $ref ${pair_id}.mapped.sorted.bam |
bcftools call --ploidy 1 -cf GQ - -Ou |
bcftools norm -f $ref - -Oz -o ${pair_id}.norm.vcf.gz