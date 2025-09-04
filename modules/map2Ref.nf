process MAP2REF {
    errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/bam", mode: 'copy', pattern: '*.bam'
	maxForks 2
	
	input:
		tuple val(pair_id), path(trim1), path(trim2)
	
	output:
    	tuple val (pair_id), path("${pair_id}.mapped.sorted.bam")
	
	script:
	"""
	#using params.ref here is an exception to normal param definition - neatest way of bringing across the index files with the ref.
	bwa mem -M -E2 -t2 ${params.ref} ${trim1} ${trim2} |
	samtools view -@2 -ShuF 2308 - |
	samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam
	"""
}