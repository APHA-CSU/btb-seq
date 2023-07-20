#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//Define variables
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

seqplate = "${params.DataDir}"
publishDir = "$params.outdir/Results_${params.DataDir}_${params.today}/"
commitId = "${workflow.commitId}"

workflow{
    /*	Collect pairs of fastq files and infer sample names
    Define the input raw sequening data files */
    Channel
        .fromFilePairs( params.reads, flat: true )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	    .set { readPairs }

    // Collect name of data folder and analysis run date
    FirstFile = file( params.reads ).first()
	    DataPath = FirstFile.getParent()
	    GetTopDir = ~/\/\*\//
	    TopDir = DataPath - GetTopDir
	    params.DataDir = TopDir.last()
	    params.today = new Date().format('ddMMMYY')

    deduplicate(readPairs)

    trim(deduplicate.out)

    map2Ref(trim.out)

    mask(map2Ref.out)

    varCall(map2Ref.out)

    readStats(read_pairs, deduplicate.out, trim.out, map2Ref.out)

    vcf2Consensus(varCall.out, mask.out)

    assignCluster(varCall.out, readStats.out)

    idNonBovis(trim.out, readStats.out)

    combineOutput(assignCluster.out, idNonBovis.out)
}
