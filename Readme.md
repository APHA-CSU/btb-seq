# **BovTB-nf**

[![APHA-CSU](https://circleci.com/gh/APHA-CSU/BovTB-nf.svg?style=svg)](https://app.circleci.com/pipelines/github/APHA-CSU)

------------

This is the updated pipeline for APHA's processing of *Mycobacterium bovis* WGS data. BovTB-nf is designed to process a batch (1 or more samples) of paired-end fastq files generated on an Illumina sequencer. It will first remove duplicate reads from the dataset (FastUniq) and then trim the unique reads based on base-call quality and the presence of adapters (Trimmomatic). Reads are then mapped to the *M. bovis* AF2122 reference genome and variants called (bwa/samtools/bcftools).

It has been built to run using nextflow, using standard bioinformatic tools for the most part. The external dependancies are:
-	FastUniq
-	Trimmomatic
-	bwa
-	samtools and bcftools
-	bedtools
-	Kraken2 (and database)
-	Bracken

## Installation

Of course Nextflow itself is a prerequisite and should be installed as described in the [Nextflow Documentation](https://www.nextflow.io/docs/latest/getstarted.html)

If you have the dependancies installed the pipeline can run by simply typing: 

	nextflow run ellisrichardj/BovTB-nf

Alternatively, clone the repository:

	git clone https://github.com/ellisrichardj/BovTB-nf.git

If required, there is simple script for installing the dependancies (helpfully called Install_dependancies.sh), which will also update the nextflow config file with their locations.

# Docker build

Alternatively, the pipeline can run in a docker ubuntu image. 

To pull the latest image (only downloads if it's not already fetched) and run the nextflow container on data:
```
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/
```

Or, build and run the image directly from source
```
docker build /PATH/TO/REPO/ -t my-bov-tb
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/ my-bov-tb
```

# Tests

Integration tests can execute in a local container built from the local `Dockerfile`
```
chmod +x tests/run_tests
./tests/run_tests
```

To add a test case, add a bash script under `tests/jobs/` and update the `JOBS` variable in `tests/run_tests`. Any command that exits with a non-zero exit code results in failure. 

# Examples

In its simplest form just run the Nextflow process from the directory containing the fastq files:

	cd /path/to/Data
	nextflow run BovTB-nf
