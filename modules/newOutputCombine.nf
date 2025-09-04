process NEWOUTPUTCOMBINE {
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}", mode: 'copy', pattern: '*.csv'
	tag "Combining outputs"
	
	input:
		path assigned
        path stats
		path qbovis
		path ncount
		val DataDir
		val user

	output:
		path('*.csv')
	
	script:
	"""
	newcombineCsv.py ${assigned} ${stats} ${qbovis} ${ncount} ${DataDir} ${workflow.commitId} ${user}
	"""
}