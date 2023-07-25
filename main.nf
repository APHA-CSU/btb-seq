#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Define variables */
params.lowmem = ""
params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"

ref = file(params.ref)
refgbk = file(params.refgbk)
rptmask = file(params.rptmask)
allsites = file(params.allsites)
stage1pat = file(params.stage1pat)
adapters = file(params.adapters)
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
	    tuple val(pair_id), file("pair_1"), file("pair_2")
	output:
	    tuple val(pair_id), file("dedup_1.fastq"), file("dedup_2.fastq")
	"""
	deduplicate.bash $pair_1 $pair_2 dedup_1.fastq dedup_2.fastq
	"""
}

process trim {
    errorStrategy 'finish'
    tag "$pair_id"
	maxForks 2
	input:
	    tuple val(pair_id), file("read_1.fastq"), file("read_2.fastq")
	output:
	    tuple val(pair_id), file("trimmed_1.fastq"), file("trimmed_2.fastq")
	"""
	trim.bash $adapters read_1.fastq read_2.fastq trimmed_1.fastq trimmed_2.fastq
	"""
}

process map2Ref {
    errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$publishDir/bam", mode: 'copy', pattern: '*.bam'
	maxForks 2
	input:
    	tuple val(pair_id), file("read_1.fastq"), file("read_2.fastq")
	output:
    	tuple val(pair_id), file("${pair_id}.bam")
	"""
	map2ref.bash $ref read_1.fastq read_2.fastq ${pair_id}.bam
	"""
}

process varCall {
	errorStrategy 'finish'
	tag "$pair_id"
	publishDir "$publishDir/vcf", mode: 'copy', pattern: '*.vcf.gz'
	maxForks 3
	input:
		tuple val(pair_id), file("mapped.bam")
	output:
		tuple val(pair_id), file("${pair_id}.vcf.gz"), file("${pair_id}.vcf.gz.csi")
	"""
	varCall.bash $ref mapped.bam ${pair_id}.vcf.gz $params.MAP_QUAL $params.BASE_QUAL $params.PLOIDY
	"""
}

workflow{
    /*	Collect pairs of fastq files and infer sample names
    Define the input raw sequening data files */
    Channel
        .fromFilePairs( params.reads, flat: true )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	    .set { readPairs }

    deduplicate(readPairs)

    trim(deduplicate.out)

    map2Ref(trim.out)

	varCall(map2Ref.out)
	
	/*mask(map2Ref.out)

    readStats(read_pairs, deduplicate.out, trim.out, map2Ref.out)

    vcf2Consensus(varCall.out, mask.out)

    assignCluster(varCall.out, readStats.out)

    idNonBovis(trim.out, readStats.out)

    combineOutput(assignCluster.out, idNonBovis.out)*/
}
