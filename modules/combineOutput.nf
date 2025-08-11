process COMBINEOUTPUT {
	tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}", mode: 'copy', pattern: '*.csv'
	
	input:
		tuple val(pair_id), path(assigned)
		tuple val(pair_id), path(qbovis)
		tuple val(pair_id), path(ncount)
		val (DataDir) 
		val (user)
	
	output:
		tuple val(pair_id), path('*.csv')
	
	script:
	"""
	combineCsv.py ${assigned} ${qbovis} ${ncount} ${DataDir} ${workflow.commitId} ${user}
	"""
}