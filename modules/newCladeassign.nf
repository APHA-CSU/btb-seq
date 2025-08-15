process NEWCLADEASSIGN {
	tag "$pair_id"
	maxForks 1
	
	input:
		tuple val(pair_id), path(consensus), path (snps)
		path (css)
	
	output:
		path("${pair_id}_cladematch.csv")

	script:
	"""
	assignClade.py ${consensus} ${css} > "${pair_id}_cladematch.csv"
	"""
}