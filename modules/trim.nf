process TRIM {
    errorStrategy 'finish'
    tag "$pair_id"
	maxForks 2
    cleanup = false
	
	input:
		tuple val(pair_id), path(uniq)
        val adapters
	
	output:
	    tuple val(pair_id), path("*trim_*.fastq")
	
	script:
	"""
	java -jar /usr/local/bin/trimmomatic.jar PE \
    -threads 2 \
    -phred33 \
    ${uniq[0]} \
    ${uniq[1]} \
    ${pair_id}_trim_R1.fastq \
    fail1.fastq \
    ${pair_id}_trim_R2.fastq \
    fail2.fastq \
    ILLUMINACLIP:${adapters}:2:30:10 \
    SLIDINGWINDOW:10:20 \
    MINLEN:36

	rm fail1.fastq fail2.fastq
	"""
}