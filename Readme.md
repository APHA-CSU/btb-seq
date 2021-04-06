# **btb-seq**

[![APHA-CSU](https://circleci.com/gh/APHA-CSU/BovTB-nf.svg?style=svg)](https://app.circleci.com/pipelines/github/APHA-CSU)

`btb-seq` is the Whole Genome Sequencing Pipeline for APHA's processing of *Mycobacterium bovis* WGS data. This pipeline is designed to process a batch (1 or more samples) of paired-end fastq files generated on an Illumina sequencer. It firsts remove duplicate reads from the dataset using `FastUniq` and then trims the unique reads based on base-call quality and the presence of adapters with `Trimmomatic`. Reads are then mapped to the *M. bovis* AF2122 reference genome and variants called (`bwa`/`samtools`/`bcftools`).

It has been built to run using [nextflow](https://www.nextflow.io/docs/latest/getstarted.html), using standard bioinformatic tools for the most part. The external dependancies are:
-	FastUniq
-	Trimmomatic
-	bwa
-	samtools and bcftools
-	bedtools
-	Kraken2 (and database)
-	Bracken


## Run pipeline in Docker

We recommend running the pipeline with our docker image. 

To pull the latest image (if it's not already fetched) and run the nextflow container on data:
```
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/
```

The `/PATH/TO/READS/` directory should contain fastq files named with the `*_{S*_R1,S*_R2}*.fastq.gz` pattern. For example, a directory with `bovis_S1_R1.fastq.gz` and `bovis_S1_R2.fastq.gz` contains a single pair of reads


### Build image from source 

You can also build and run the image directly from source
```
docker build /PATH/TO/REPO/ -t my-bov-tb
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/ my-bov-tb
```


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

