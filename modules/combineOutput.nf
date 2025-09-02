process COMBINEOUTPUT {
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}", mode: 'copy', pattern: '*.csv'
	tag "Combining outputs"
	
	input:
		path assigned
		path qbovis
		path ncount
		path DataDir
		val user

	output:
		path('*.csv')
	
	script:
	"""
	combineCsv.py ${assigned} ${qbovis} ${ncount} ${DataDir} ${workflow.commitId} ${user}
	"""
}