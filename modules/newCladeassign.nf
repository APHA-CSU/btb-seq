process NEWCLADEASSIGN {
	tag "$pair_id"
	maxForks 1
	
	input:
		tuple val(pair_id), path(consensus), path (snps)
		path params.csstable
	
	output:
		path("${pair_id}_cladematch.csv")

	script:
	"""
	assignClade.py ${consensus} ${params.csstable} ${pair_id}_cladematch.csv
	"""
}