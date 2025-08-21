process DEDUPLICATE {
    errorStrategy 'finish'
    tag "$pair_id"
	maxForks 2
	
	input:
        tuple val(pair_id), path(read1), path(read2)
	
	output:
        tuple val(pair_id), path("${pair_id}_uniq_R1.fastq"), path("${pair_id}_uniq_R2.fastq")
	
	script:
	"""
	set -eo pipefail

	# Unzip .fastq.gz
	gunzip -c ${read1} > ${pair_id}_R1.fastq 
	gunzip -c ${read2} > ${pair_id}_R2.fastq

	# Input List
    echo -e "${pair_id}_R1.fastq\\n${pair_id}_R2.fastq" > fqin.lst

	# Extract unique reads
	fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
	"""
}