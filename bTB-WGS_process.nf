#!/usr/bin/env nextflow

/*	This is APHA's nextflow pipeline to process Illumina paired-end data from Mycobacterium bovis isolates.  
*	It will first deduplicate the reads using fastuniq, trim them using trimmomatic and then map to the reference genome.
*	Variant positions wrt the reference are determined, togther with data on the number of reads mapping and the depth of 
*	coverage.  Using a panel of predetermined cluster specific SNPs it will also infer cluster membership.
*/

/* Default parameters */
params.lowmem = ""
params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"

ref = file(params.ref)
refgbk = file(params.refgbk)
rptmask = file(params.rptmask)
allsites = file(params.allsites)
stage1pat = file(params.stage1pat)
stage2pat = file(params.stage2pat)
adapters = file(params.adapters)
discrimPos = file(params.discrimPos)

pypath = file(params.pypath)
kraken2db = file(params.kraken2db)

/*	Collect pairs of fastq files and infer sample names
Define the input raw sequening data files */
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs } 
	read_pairs.into { read_pairs; raw_reads }

publishDir = "$params.outdir/"

/* remove duplicates from raw data
This process removes potential duplicate data (sequencing and optical replcaites from the raw data set */
process Deduplicate {
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	tuple pair_id, pair_1, pair_2 from read_pairs

	output:
	tuple pair_id, file("dedup_1.fastq"), file("dedup_2.fastq") into dedup_read_pairs, uniq_reads

	"""
	deduplicate.bash $pair_1 $pair_2 dedup_1.fastq dedup_2.fastq
	"""
}	

/* trim adapters and low quality bases from fastq data
Removes the adapters which are added during the lab processing and and any low quality data */
process Trim {
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	tuple pair_id, file("read_1.fastq"), file("read_2.fastq") from dedup_read_pairs

	output:
	tuple pair_id, file("trimmed_1.fastq"), file("trimmed_2.fastq") into trim_read_pairs, trim_read_pairs2, trim_reads
	
	"""
	trim.bash $adapters read_1.fastq read_2.fastq trimmed_1.fastq trimmed_2.fastq
	"""
}

/* map to reference sequence
Aligns the individiual sequence reads to the reference genome */
process Map2Ref {
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$publishDir/bam", mode: 'copy', pattern: '*.bam'

	maxForks 2

	input:
	tuple pair_id, file("read_1.fastq"), file("read_2.fastq") from trim_read_pairs

	output:
	tuple pair_id, file("${pair_id}.bam") into mapped_bam, bam4stats, bam4mask

	"""
	map2Ref.bash $ref read_1.fastq read_2.fastq ${pair_id}.bam
	"""
}

/* Variant calling
Determines where the sample differs from the reference genome */
process VarCall {
	errorStrategy 'finish'
	tag "$pair_id"

	publishDir "$publishDir/vcf", mode: 'copy', pattern: '*.vcf.gz'

	maxForks 3

	input:
	tuple pair_id, file("mapped.bam") from mapped_bam

	output:
	tuple pair_id, file("${pair_id}.vcf.gz"), file("${pair_id}.vcf.gz.csi") into vcf, vcf2,	vcf4mask

	"""
	varCall.bash $ref mapped.bam ${pair_id}.vcf.gz
	"""
}

vcf4mask
	.join(bam4mask)
	.set { mask_input }

/* Masking known repeats regions and sites with zero coverage
Ensure that consensus only includes regions of the genome where there is high confidence */
process Mask {
	errorStrategy 'finish'
	tag "$pair_id"

	maxForks 2

	input:
	tuple pair_id, file("called.vcf"), file("called.vcf.csi"), file("mapped.bam") from mask_input

	output:
	tuple pair_id, file("mask.bed"), file("nonmasked-regions.bed") into maskbed

	"""
	mask.bash $rptmask called.vcf mask.bed nonmasked-regions.bed mapped.bam $allsites $params.MIN_READ_DEPTH $params.MIN_ALLELE_FREQUENCY_ALT $params.MIN_ALLELE_FREQUENCY_REF
	"""
}

// Combine input for consensus calling
// Joins relevant files for the same sample as input to the next process

maskbed
	.join(vcf2)
	.set { vcf_bed }

/* Consensus calling and output snp table for snippy compatability */
process VCF2Consensus {
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$publishDir/consensus/", mode: 'copy', pattern: '*.fas'
	publishDir "$publishDir/snpTables/", mode: 'copy', pattern: '*.tab'
	publishDir "$publishDir/filteredBcf/", mode: 'copy', pattern: '*.bcf'
	publishDir "$publishDir/filteredBcf/", mode: 'copy', pattern: '*.csi'

	maxForks 2

	input:
	tuple pair_id, file("mask.bed"), file("nonmasked-regions.bed"), file("variant.vcf.gz"),	file("variant.vcf.gz.csi") from vcf_bed

	output:
	tuple pair_id, file("${pair_id}_consensus.fas") into consensus
	tuple pair_id, file("${pair_id}_snps.tab") into snpstab
	tuple pair_id, file("${pair_id}_filtered.bcf"), file("${pair_id}_filtered.bcf.csi") into _

	"""
	vcf2Consensus.bash $ref mask.bed nonmasked-regions.bed variant.vcf.gz ${pair_id}_consensus.fas ${pair_id}_snps.tab ${pair_id}_filtered.bcf $params.MIN_ALLELE_FREQUENCY_ALT
	"""
}

//	Combine data for generating per sample statistics
// Joins relevant files for the same sample as input to the next process

raw_reads
	.join(uniq_reads)
	.set { raw_uniq }

trim_reads
	.join(bam4stats)
	.set { trim_bam }

raw_uniq
	.join(trim_bam)
	.set { input4stats }

/* Generation of data quality and mapping statistics
Calculate number of raw reads, unique reads, trimmed reads, proportion aligned to reference genome */

process ReadStats{
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	set pair_id, file("${pair_id}_*_R1_*.fastq.gz"), file("${pair_id}_*_R2_*.fastq.gz"), file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq"), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq"), file("${pair_id}.mapped.sorted.bam") from input4stats

	output:
	set pair_id, file("${pair_id}_stats.csv") into stats
	set pair_id, file('outcome.txt') into Outcome

    """
    readStats.bash "$pair_id"
    """
}

//	Combine inputs to assign cluster for each sample
// Joins relevant files for the same sample as input to the next process

vcf
	.join(stats)
	.set { input4Assign }

/* Assigns cluster by matching patterns of cluster specific SNPs
Compares SNPs identified in vcf file to lists in reference table */

process AssignClusterCSS{
	errorStrategy 'ignore'
    tag "$pair_id"
	

	maxForks 1

	input:
	set pair_id, file("${pair_id}.vcf.gz"), file("${pair_id}.vcf.gz.csi"), file("${pair_id}_stats.csv") from input4Assign

	output:
	file("${pair_id}_stage1.csv") into AssignCluster

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

// Collect data to ID any non-M.bovis samples
// Joins relevant files for the same sample as input to the next process

Outcome
	.join(trim_read_pairs2)
	.set { IDdata }

/* Identify any non-M.bovis samples using kraken
Samples with flag != 'Pass' are processed to detemine which microbe(s) are present 
Bracken parses the output which is then sorted to generate a top 20 list of species
Presence / absence of M.bovis is also determined by parsing Bracken output */

process IDnonbovis{
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$publishDir/NonBovID", mode: 'copy', pattern: '*.tab'

	maxForks 1

	input:
	set pair_id, file('outcome.txt'), file("trimmed_1.fastq"), file("trimmed_2.fastq") from IDdata

	output:
	set pair_id, file("${pair_id}_*_brackensort.tab"), file("${pair_id}_*_kraken2.tab")  optional true into IDnonbovis
	file("${pair_id}_bovis.csv") optional true into QueryBovis

	"""
	idNonBovis.bash $pair_id $kraken2db $params.lowmem
	"""
}

/* Combine all cluster assignment data into a single results file */

AssignCluster
	.collectFile( name: "${publishDir}/AssignedWgsCluster.csv", sort: true, storeDir: "$publishDir/", keepHeader: true )
	.set {Assigned}

QueryBovis
	.collectFile( name: "${publishDir}/BovPos.csv", sort: true, storeDir: "$publishDir", keepHeader: true )
	.set {Qbovis}

process CombineOutput {
	publishDir "${publishDir}/", mode: 'copy', pattern: '*.csv'
	
	input:
	file('Assigned.csv') from Assigned
	file('Qbovis.csv') from Qbovis

	output:
	file '*.csv' into FinalOut

	"""
	combineCsv.py Assigned.csv Qbovis.csv ./
	"""
}

workflow.onComplete {
		log.info "Completed sucessfully:	$workflow.success"
		log.info "Nextflow Version:	$workflow.nextflow.version"
		log.info "Duration:		$workflow.duration"
		log.info "Output Directory:	${publishDir}/"
}
