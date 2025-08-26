process NEWOUTPUTCOMBINE {
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}", mode: 'copy', pattern: '*.csv'
	tag "Combining outputs"
	
	input:
		path(assigned)
        path(stats)
		path(qbovis)
		path(ncount)

	output:
		path('*.csv')
	
	script:
	"""
	newcombineCsv.py ${assigned} ${stats} ${qbovis} ${ncount} ${params.DataDir} ${workflow.commitId} ${params.user}
	"""
}