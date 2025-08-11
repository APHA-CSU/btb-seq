process MAP2REF {
    errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/bam", mode: 'copy', pattern: '*.bam'
	maxForks 2
	cleanup = false
	
	input:
		tuple val(pair_id), path(trimmed)
        val ref
	
	output:
    	tuple val (pair_id), path ("*.mapped.sorted.bam")
	
	script:
	"""
	bwa mem -M -E 2 -t2 ${ref} ${trimmed[0]} ${trimmed[1]} |
	samtools view -@2 -ShuF 2308 - |
	samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam
	"""
}