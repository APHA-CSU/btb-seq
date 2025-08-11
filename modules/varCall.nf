process VARCALL {
	errorStrategy 'finish'
	tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/vcf", mode: 'copy', pattern: '*.vcf.gz'
	maxForks 3
	cleanup = false

	input:
		tuple val(pair_id), path(bam)
        path (ref)
        val (quality)
        val (base_qual)
        val (ploidy)
	
	output:
		tuple val (pair_id), path ("*vcf.gz"), path ("*csi")
	
	script:
	"""
	samtools index ${bam}
	bcftools mpileup -q ${quality} -Q ${base_qual} -a INFO/AD -Ou -f ${ref} ${bam} |
    bcftools call --ploidy ${ploidy} -mf GQ - -Ou |
    bcftools norm -f ${ref} - -Oz -o ${pair_id}.vcf.gz
	bcftools index ${pair_id}.vcf.gz
	"""
}