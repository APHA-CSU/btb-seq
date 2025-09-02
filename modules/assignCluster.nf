process ASSIGNCLUSTER {
	errorStrategy 'ignore'
    tag "$pair_id"
	maxForks 1
	
	input:
		tuple val(pair_id), path(vcf), path(csi), path(stats)
		path discrimPos 
		path stage1pat
		path ref
		val min_mean_cov 
		val min_cov_snp
		val alt_prop_snp
		val min_qual_snp
		val min_qual_nonsnp
		path pypath
	
	output:
		path("${pair_id}_stage1.csv")
	
	script:
	"""
	set -eo pipefail

	bcftools view -R ${discrimPos} -O v -o ${pair_id}.discrimPos.vcf ${vcf}
	python3 ${pypath}/Stage1-test.py ${stats} ${stage1pat} ${ref} test ${min_mean_cov} ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} ${min_qual_nonsnp} ${pair_id}.discrimPos.vcf
	mv _stage1.csv ${pair_id}_stage1.csv
	"""
}