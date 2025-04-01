#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Define variables */
params.lowmem = ""
params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"

//ref = file(params.ref)
//refgbk = file(params.refgbk)
//rptmask = file(params.rptmask)
allsites = file(params.allsites)
stage1pat = file(params.stage1pat)
//adapters = file(params.adapters)
discrimPos = file(params.discrimPos)

pypath = file(params.pypath)
kraken2db = file(params.kraken2db)

//Collect name of data folder and analysis run date
    FirstFile = file( params.reads ).first()
	DataPath = FirstFile.getParent()
	GetTopDir = ~/\/\*\//
	TopDir = DataPath - GetTopDir
	params.DataDir = TopDir.last()
	params.today = new Date().format('ddMMMYY')

	seqplate = "${params.DataDir}"
	publishDir = "$params.outdir/Results_${params.DataDir}_${params.today}/"
	commitId = "${workflow.commitId}"

process deduplicate {
    errorStrategy 'finish'
    tag "$pair_id"
	maxForks 2
	input:
	    tuple val(pair_id), path("pair_1"), path("pair_2")
	output:
	    tuple val(pair_id), path("dedup_1.fastq"), path("dedup_2.fastq")
	"""
	deduplicate.bash $pair_1 $pair_2 dedup_1.fastq dedup_2.fastq
	"""
}

process trim {
    errorStrategy 'finish'
    tag "$pair_id"
	maxForks 2
	input:
	    tuple val(pair_id), path("read_1.fastq"), path("read_2.fastq")
		path("adapters.fas")
	output:
	    tuple val(pair_id), path("trimmed_1.fastq"), path("trimmed_2.fastq")
	"""
	trim.bash adapters.fas read_1.fastq read_2.fastq trimmed_1.fastq trimmed_2.fastq
	"""
}

process map2Ref {
    errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$publishDir/bam", mode: 'copy', pattern: '*.bam'
	maxForks 2
	input:
    	tuple val(pair_id), path("read_1.fastq"), path("read_2.fastq")
		path("ref.fas")
	output:
    	tuple val(pair_id), path("${pair_id}.bam")
	script:
	"""
	map2Ref.bash ref.fas read_1.fastq read_2.fastq ${pair_id}.bam
	"""
}

process varCall {
	errorStrategy 'finish'
	tag "$pair_id"
	publishDir "$publishDir/vcf", mode: 'copy', pattern: '*.vcf.gz'
	maxForks 3
	input:
		tuple val(pair_id), path("mapped.bam")
		path("ref.fas")
	output:
		tuple val(pair_id), path("${pair_id}.vcf.gz"), path("${pair_id}.vcf.gz.csi")
	script:
	"""
	varCall.bash ref.fas mapped.bam ${pair_id}.vcf.gz $params.MAP_QUAL $params.BASE_QUAL $params.PLOIDY
	"""
}

process mask {
	errorStrategy 'finish'
	tag "$pair_id"
	maxForks 2
	input:
		tuple val(pair_id), path("called.vcf"), path("called.vcf.csi"), path("mapped.bam")
		path("rptmask.bed")
	output:
		tuple val(pair_id), path("mask.bed"), path("nonmasked-regions.bed")
	script:
	"""
	mask.bash rptmask.bed called.vcf mask.bed nonmasked-regions.bed mapped.bam $allsites $params.MIN_READ_DEPTH $params.MIN_ALLELE_FREQUENCY_ALT $params.MIN_ALLELE_FREQUENCY_REF
	"""
}

process readStats {
	errorStrategy 'finish'
    tag "$pair_id"
	maxForks 2
	input:
		tuple val(pair_id), path("${pair_id}_*_R1_*.fastq.gz"), path("${pair_id}_*_R2_*.fastq.gz"), 
		path("${pair_id}_uniq_R1.fastq"), path("${pair_id}_uniq_R2.fastq"),
		path("${pair_id}_trim_R1.fastq"), path("${pair_id}_trim_R2.fastq"),
		path("${pair_id}.mapped.sorted.bam")
	output:
		tuple val(pair_id), path("${pair_id}_stats.csv"), emit: stats
		tuple val(pair_id), path('outcome.txt'), emit: outcome
    """
    readStats.bash "$pair_id"
    """
}

process vcf2Consensus {
	errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$publishDir/consensus/", mode: 'copy', pattern: '*.fas'
	publishDir "$publishDir/snpTables/", mode: 'copy', pattern: '*.tab'
	maxForks 2
	input:
		tuple val(pair_id), path("mask.bed"), path("nonmasked-regions.bed"),
		path("variant.vcf.gz"),	path("variant.vcf.gz.csi")
		path("ref.fas")
	output:
		tuple val(pair_id), path("${pair_id}_consensus.fas"), path("${pair_id}_snps.tab"), emit: consensus
		path("${pair_id}_ncount.csv"), emit: nCount
	script:
	"""
	vcf2Consensus.bash ref.fas \
		mask.bed \
		nonmasked-regions.bed \
		variant.vcf.gz \
		${pair_id}_consensus.fas \
		${pair_id}_snps.tab \
		${pair_id}_filtered.bcf \
		$params.MIN_ALLELE_FREQUENCY_ALT \
		$pair_id \
		$publishDir
	"""
}

process assignCluster {
	errorStrategy 'ignore'
    tag "$pair_id"
	maxForks 1
	input:
		tuple val(pair_id), path("${pair_id}.vcf.gz"), path("${pair_id}.vcf.gz.csi"),
		path("${pair_id}_stats.csv")
	output:
		path("${pair_id}_stage1.csv")
	"""
	assignClusterCss.bash $pair_id \
		$discrimPos \
		$stage1pat \
		$min_mean_cov \
		$min_cov_snp \
		$alt_prop_snp \
		$min_qual_snp \
		$min_qual_nonsnp \
		$pypath \
		$ref
	"""
}

process idNonBovis {
	errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/NonBovID", mode: 'copy', pattern: '*.tab'
	maxForks 1
	input:
		tuple val(pair_id), path("outcome.txt"), path("trimmed_1.fastq"), path("trimmed_2.fastq")
	output:
		path("${pair_id}_bovis.csv"), emit: queryBovis
		tuple val(pair_id), path("${pair_id}_*_brackensort.tab"), path("${pair_id}_*_kraken2.tab"), optional: true, emit: krakenOut
	"""
	idNonBovis.bash $pair_id $kraken2db $params.lowmem
	"""
}

process combineOutput {
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}", mode: 'copy', pattern: '*.csv'
	input:
		path('assigned_csv')
		path('qbovis_csv')
		path('ncount_csv')
	output:
		path('*.csv')
	"""
	combineCsv.py assigned_csv qbovis_csv ncount_csv $seqplate $commitId
	"""
}

workflow{
    /*	Collect pairs of fastq files and infer sample names
    Define the input raw sequening data files */
    Channel
        .fromFilePairs( params.reads, flat: true )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	    .set { readPairs }

	Channel
		.fromPath( params.adapters )
		.set { adapters }

	Channel
		.fromPath( params.ref )
		.set { ref }

	Channel
		.fromPath( params.rptmask)
		.set { rptmask }

	deduplicate(readPairs)

	trim(deduplicate.out, adapters)

	map2Ref(trim.out, ref)

	varCall(map2Ref.out, ref)
	
	varCall.out
		.join( map2Ref.out )
		.set {vcf_bam}

	mask(vcf_bam, rptmask)

	readPairs
		.join( deduplicate.out )
		.join( trim.out )
		.join( map2Ref.out )
		.set {reads_mapped}

	readStats(reads_mapped)

	mask.out
		.join( varCall.out )
		.set {mask_vcf}

	vcf2Consensus(mask_vcf, ref)

	varCall.out
		.join( readStats.out.stats )
		.set {vcf_stats}

	assignCluster(vcf_stats)

	readStats.out.outcome
		.join( trim.out )
		.set {outcome_reads}

	idNonBovis(outcome_reads)

	assignCluster.out
		.collectFile( name: "${params.DataDir}_AssignedWGSCluster_${params.today}.csv", sort: true, keepHeader: true )
		.set {assigned}

	idNonBovis.out.queryBovis
		.collectFile( name: "${params.DataDir}_BovPos_${params.today}.csv", sort: true, keepHeader: true )
		.set {qbovis}

	vcf2Consensus.out.nCount
		.collectFile( name: "${params.DataDir}_Ncount_${params.today}.csv", sort: true, keepHeader: true)
		.set {consensusQual}

	combineOutput(assigned, qbovis, consensusQual)
}
