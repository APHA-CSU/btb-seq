#!/usr/bin/env nextflow

/*	This is APHA's nextflow pipeline to process Illumina paired-end data from Mycobacterium bovis isolates.  
*	It will first deduplicate the reads using fastuniq, trim them using trimmomatic and then map to the reference genome.
*	Variant positions wrt the reference are determined, togther with data on the number of reads mapping and the depth of 
*	coverage.  Using a panel of predetermined cluster specific SNPs it will also infer cluster membership.
*
*	written by ellisrichardj, based on pipeline developed by Javier Nunez
*
*	Version 0.1.0	31/07/18	Initial version
*	Version 0.2.0	04/08/18	Added SNP filtering and annotation
*	Version 0.3.0	06/08/18	Generate summary of samples in batch
*	Version 0.4.0	16/08/18	Minor improvements
*	Version 0.5.0	14/09/18	Infer genotypes using genotype-specific SNPs (GSS)
*	Version 0.5.1	21/09/18	Fixed bug that prevented GSS from running correctly
*	Version 0.5.2	01/10/18	Mark shorter split hits as secondary in bam file (-M) and change sam flag filter to 3844
*	Version 0.5.3	15/10/18	Increase min mapQ for mpileup to 60 to ensure unique reads only; add $dependpath variable
*	Version 0.6.0	15/11/18	Assigns clusters (newly defined) in place of inferring historical genotype.
*	Version 0.6.1	15/11/18	Fixed bug which caused sample names to be inconsistently transferred between processes
*	Version 0.6.2	24/11/18	Fixed bug to allow cluster assignment to be collected in a single file
*	Version 0.6.3	26/11/18	Removed 'set' in output declaration as it caused nextflow warning
*	Version 0.7.0	26/11/18	Add process to output phylogenetic tree
*	Version 0.7.1	11/12/18	Used join to ensure inputs are properly linked
*	Version 0.7.2	18/12/18	Changed samtools filter to remove unmapped reads and reads that aligned more than once
*	Version 0.7.3	22/01/19	Changed link to adapters file for trimming
*	Version 0.7.4	22/02/19	Changed calculations of mapping stats
*	Version 0.7.5	08/03/19	Further correction of mapping stats and addition of success flags
*	Version 0.8.0	12/03/19	Add kraken to ID samples that fail cluster assignment
*	Version 0.8.1	13/03/19	Update to kraken2
*	Version 0.8.2	14/03/19	Define location of kraken2 database as a nextflow parameter
*	Version 0.8.3	14/03/19	Add option to reduce memory use by kraken2 if required
*	Verison 0.8.4	15/03/19	Correct output location of kraken2 tables
*	Version 0.8.5	25/03/19	Output bam, vcf and consensus to Results directory
*	Version 0.8.6	26/03/19	Add normalization and filtering of indels in vcf before consensus calling
*	Version 0.8.7	29/04/19	Remove intermediary fastq files, rebalance processes and remove redundant process
*	Version 0.8.8	03/05/19	More rebalancing and removing redundant output
*	Version 0.9.0	10/09/19	Filter and mask vcf for consensus calling
*	Version 0.9.1	19/09/19	Remove SNP filtering and annotation process as no longer required
*	Version 0.9.2	20/09/19	Exclude indels from consensus calling step
*	Version 0.9.3   11/10/19	Move ReadStats to standalone shell script
*	Version 0.9.4	17/10/19	Add Data source and datestamp to output files and directories
*	Version 0.9.5	20/12/19	Update to Python3 and lastest samtools/bcftools (1.10)
*	Version 0.9.6	30/03/20	Add Bracken for parsing kraken2 output
*	Version 0.9.7	21/04/20	Determine presence of M.bovis in low quality samples
*	Version 0.9.8	30/04/20	Output snp table in tsv format for each sample
*/

process test{

	publishDir "$params.outdir", mode: 'copy'

	container "aaronsfishman/bov-tb:batch-no-deduplicate"

	scratch '/BovTB-nf/nxf-scratch-dir'

	output:
		file "scratch_test.txt"

	script:
	"""
	touch scratch_test.txt
	echo \"PWD `pwd` \" >> scratch_test.txt
	"""
}
