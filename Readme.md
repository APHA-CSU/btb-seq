# **btb-seq**

[![APHA-CSU](https://circleci.com/gh/APHA-CSU/btb-seq.svg?style=svg)](https://app.circleci.com/pipelines/github/APHA-CSU)

`btb-seq` is the Whole Genome Sequencing Pipeline for APHA's processing of raw *Mycobacterium bovis* WGS data. This pipeline is designed to process a batch (1 or more samples) of paired-end `fastq.gz` files generated on an Illumina sequencer using [nextflow](https://www.nextflow.io/docs/latest/getstarted.html). 

## Installation

To install the software in ubuntu, run
```
bash install.bash
```

This script installs the following dependancies and adds symlinks to the `$PATH`: 
-	`nextflow`
-	`FastUniq`
-	`Trimmomatic`
-	`bwa`
-	`samtools` and `bcftools`
-	`bedtools`
-	`Kraken2` (and database)
-	`Bracken`

## Running the pipeline

To run the pipeline on a batch a samples, a directory containing raw `.fastq.gz` files is required. Each read-pair sample is represented by a pair of files named `*_R1.fastq.gz` and `*_R2.fastq.gz`. For example, to batch two samples named `bovis_a` and `bovis_b`, a directory containing `bovis_a_R1.fastq.gz`, `bovis_a_R2.fastq.gz`,  `bovis_b_R1.fastq.gz` and `bovis_b_R2.fastq.gz`, needs to be defined.

Pipeline output is stored in a results directory that contains
- Consensus `fasta`
- `.bam` files
- Kraken metagenomics `.tab`
- A summary csv file
- ... 

### Run from terminal

To run a batch from the terminal
```
./btb-seq /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/
```

### Run from docker

To run the pipeline from ther terminal call
```
TODO
```
This will pull the latest image (if it's not already fetched) and run the nextflow container on data.

### Run from docker

```
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/
```


### Build image from source 

You can also build your own experimental docker fom source
```
docker build ./docker/ -t my-bov-tb
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/ my-bov-tb
```

## How it works

Duplicate reads are removed from the dataset using `FastUniq` and then trims the unique reads based on base-call quality and the presence of adapters with `Trimmomatic`. Reads are then mapped to the *M. bovis* AF2122 reference genome and variants called (`bwa`/`samtools`/`bcftools`).


## Validation

The pipeline is validated against real-world biological samples sequenced with Illumina NextSeq machines at APHA. The test code for the validation tests is stored under `tests/jobs/`. A summary of each test is described below


### Quality Test

The quality test ensures that low quality reads (<20) are not considered for variant calling and genotyping. This is performed by setting uniform quality values to a real-world *M. bovis* sample and asserting output. Low quality bases are removed from the sequence using `Trimmomatic`, which uses a sliding window that deletes reads when the average base quality drops below 20. A table of expected results is shown below.

| Base Quality | Outcome | flag | group |
| ------------- | ------------- | ------------- | ------------- | 
| 19   | CheckRequired | LowCoverage | NA |
| 20   | Pass | BritishbTB | B6-16 |

### Limit of Detection (LoD)

The limit of detection test ensures mixtures of M. Avium and M. Bovis at varying proportions give the correct Outcome. This is performed by taking random reads from reference samples of M. Bovis and M. Avium.


| M. Bovis (%) | M. Avium (%) | Outcome |
| ------------- | ------------- | ------------- | 
| 100%   | 0% | Pass | 
| 65%   | 35% | BritishbTB | 
| 60%   | 40% | CheckRequired | 
| 0%   | 100% | Comtaminated | 

