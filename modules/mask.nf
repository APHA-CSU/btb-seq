process MASK {
	errorStrategy 'finish'
	tag "$pair_id"
	maxForks 2
	
	input:
        tuple val(pair_id), path(vcf), path(csi)
        tuple val(pair_id), path(bam)
        path (rptmask_bed)
        path (allsites_bed)
        val (min_read_depth)
        val (min_freq_alt) 
        val (min_freq_ref)
	
	output:
        tuple val(pair_id), path("${pair_id}.mask.bed"), path("${pair_id}.nonmasked-regions.bed")
	
	script:
	"""
	# Construct a mask: 
	# mask regions which don't have {sufficient evidence for alt AND sufficient evidence for the REF}
	bcftools filter -i "(ALT!='.' && INFO/AD[1] < ${min_read_depth} && INFO/AD[0]/(INFO/AD[0]+INFO/AD[1]) <= ${min_freq_ref}) ||
    (ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) < ${min_freq_alt} && INFO/AD[0]/(INFO/AD[0]+INFO/AD[1]) <= ${min_freq_ref})" ${vcf} -ov -o ${pair_id}.quality-mask.vcf
	bedtools merge -i ${pair_id}.quality-mask.vcf > ${pair_id}.quality-mask.bed

	# mask regions where there is zero coverage
	bedtools genomecov -bga -ibam ${bam} | (grep -w "0\$" || true) | cut -f -3 > ${pair_id}.zerocov.bed

	# Merge with exisiting known repeat regions
	cat ${pair_id}.quality-mask.bed ${pair_id}.zerocov.bed ${rptmask_bed} | 
	sort -k1,1 -k2,2n |
	bedtools merge > ${pair_id}.mask.bed

	# Make bedfile of sites to keep
	bedtools subtract -a ${allsites_bed} -b ${pair_id}.mask.bed > ${pair_id}.nonmasked-regions.bed

	# Cleanup
	rm ${pair_id}.quality-mask.vcf ${pair_id}.quality-mask.bed ${pair_id}.zerocov.bed
	"""
}