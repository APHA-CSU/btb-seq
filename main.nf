#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* Define variables */
params.lowmem = ""
params.reads = "${env('PWD')}/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "${env('PWD')}"
params.help = false

/* Define run information */
params.DataDir = params.reads.tokenize('/')[-2]
params.today = new Date().format('ddMMMYY')
params.user = "UnknownUser"

/* location of reference information */
params.ref = "$projectDir/references/Mycbovis-2122-97_LT708304.fas"
params.refgbk = "$projectDir/references/Mycbovis-2122-97_LT708304.gb"
params.rptmask = "$projectDir/references/DataDrivenMerge20.bed"
params.allsites = "$projectDir/references/All-sites.bed"
params.adapters = "$projectDir/references/adapter.fasta"
params.discrimPos = "$projectDir/references/DiscrimPos.tsv"
params.csstable = "$projectDir/references/CSSwithref.csv"

/* params for original assignCluster process */
params.stage1pat = "$projectDir/references/Stage1_patterns/"
params.pypath = "$projectDir/pyscripts/"

/* quality thresholds */
params.min_mean_cov = 8
params.min_cov_snp = 5
params.alt_prop_snp = 0.2
params.min_qual_snp = 150
params.min_qual_nonsnp = 0

/* quality thresholds (variant calling) */
params.MAP_QUAL = 0
params.BASE_QUAL = 10
params.PLOIDY = "1"

/* quality thresholds (snp filtering) */
params.MIN_READ_DEPTH = 8
params.MIN_ALLELE_FREQUENCY_ALT = 0.8
params.MIN_ALLELE_FREQUENCY_REF = 0.7

/* location of dependancies */
params.kraken2db = "/opt/Kraken2/db/minikraken2_v1_8GB/"

/* Print help message if --help is passed */
if (params.help) {
    println """
    Usage: nextflow run main.nf [options]

    Options:
      --reads                  Path to input FASTQ files (default: ${params.reads})
      --outdir                 Output directory (default: ${params.outdir})
      --lowmem                 Enable low memory mode for Kraken2 (default: ${params.lowmem})
      --kraken2db              Path to Kraken2 database (default: ${params.kraken2db})
      --ref                    Path to reference genome (default: ${params.ref})
      --refgbk                 Path to reference GenBank file (default: ${params.refgbk})
      --rptmask                Path to repeat mask file (default: ${params.rptmask})
      --allsites               Path to all-sites BED file (default: ${params.allsites})
      --adapters               Path to adapter sequences file (default: ${params.adapters})
      --discrimPos             Path to discrimination positions file (default: ${params.discrimPos})
      --csstable               Path to CSS table file (default: ${params.csstable})
      --stage1pat              Path to Stage1 patterns directory (default: ${params.stage1pat})
      --pypath                 Path to Python scripts directory (default: ${params.pypath})
      --min_mean_cov           Minimum mean coverage (default: ${params.min_mean_cov})
      --min_cov_snp            Minimum coverage for SNPs (default: ${params.min_cov_snp})
      --alt_prop_snp           Minimum alternate proportion for SNPs (default: ${params.alt_prop_snp})
      --min_qual_snp           Minimum quality for SNPs (default: ${params.min_qual_snp})
      --min_qual_nonsnp        Minimum quality for non-SNPs (default: ${params.min_qual_nonsnp})
      --MAP_QUAL               Minimum mapping quality (default: ${params.MAP_QUAL})
      --BASE_QUAL              Minimum base quality (default: ${params.BASE_QUAL})
      --PLOIDY                 Ploidy value for variant calling (default: ${params.PLOIDY})
      --MIN_READ_DEPTH         Minimum read depth for SNP filtering (default: ${params.MIN_READ_DEPTH})
      --MIN_ALLELE_FREQUENCY_ALT Minimum alternate allele frequency (default: ${params.MIN_ALLELE_FREQUENCY_ALT})
      --MIN_ALLELE_FREQUENCY_REF Minimum reference allele frequency (default: ${params.MIN_ALLELE_FREQUENCY_REF})
      --help                   Print this help message and exit

    Example:
      nextflow run main.nf --reads '/path/to/data/*_{S*_R1,S*_R2}*.fastq.gz' -with-docker aphacsubot/btb-seq:master'
    """
    exit 0
}

/* set modules */
include { DEDUPLICATE } from './modules/deduplicate'
include { TRIM } from './modules/trim'
include { MAP2REF } from './modules/map2Ref'
include { VARCALL } from './modules/varCall'
include { MASK } from './modules/mask'
include { VCF2CONSENSUS } from './modules/vcf2Consensus'
include { ASSIGNCLUSTER  } from './modules/assignCluster'
include { NEWCLADEASSIGN } from './modules/newCladeassign'
include { READSTATS } from './modules/readStats'
include { IDNONBOVIS } from './modules/idNonBovis'
include { COMBINEOUTPUT } from './modules/combineOutput'

/* define workflow */
workflow {
	ch_reads = Channel.fromFilePairs(params.reads, flat: true).ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	deduplicate_ch = DEDUPLICATE(ch_reads)
  trim_ch = TRIM(deduplicate_ch, params.adapters)
	map2Ref_ch = MAP2REF(trim_ch, params.ref)
	varCall_ch = VARCALL(map2Ref_ch, params.ref, params.MAP_QUAL, params.BASE_QUAL, params.PLOIDY)
	mask_ch = MASK(varCall_ch, map2Ref_ch, params.rptmask, params.allsites, params.MIN_READ_DEPTH, params.MIN_ALLELE_FREQUENCY_ALT, params.MIN_ALLELE_FREQUENCY_REF)
	vcf2Consensus_ch = VCF2CONSENSUS(params.ref, mask_ch.join(varCall_ch), params.MIN_ALLELE_FREQUENCY_ALT, params.outdir, params.today)
	readStats_ch = READSTATS(ch_reads.join(deduplicate_ch).join(trim_ch).join(map2Ref_ch))
	assignCluster_ch = ASSIGNCLUSTER(varCall_ch.join(readStats_ch.stats), params.discrimPos, params.stage1pat, params.ref, params.min_mean_cov, params.min_cov_snp, params.alt_prop_snp, params.min_qual_snp, params.min_qual_nonsnp, params.pypath)
	newcladeassign_ch = NEWCLADEASSIGN(vcf2Consensus_ch.consensus, params.csstable)
  idNonBovis_ch = IDNONBOVIS(readStats_ch.outcome.join(trim_ch), params.kraken2db, params.lowmem)
  assignCluster_ch.collectFile(name: "${params.DataDir}_AssignedWGSCluster_${params.today}.csv", sort: true, keepHeader: true).set{assigned}
  newcladeassign_ch.collectFile(name: "${params.DataDir}_AssignedClade_${params.today}.csv", keepHeader: true, storeDir: "${params.outdir}/Results_${params.DataDir}_${params.today}").set{newclade}
  idNonBovis_ch.queryBovis.collectFile(name: "${params.DataDir}_BovPos_${params.today}.csv", sort: true, keepHeader: true).set {qbovis}
  vcf2Consensus_ch.nCount.collectFile(name: "${params.DataDir}_Ncount_${params.today}.csv", sort: true, keepHeader: true).set {consensusQual}
  combineoutput_ch = COMBINEOUTPUT(assigned, qbovis, consensusQual, params.DataDir, params.user)
}