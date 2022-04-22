# **btb-seq**

<img src="https://user-images.githubusercontent.com/6979169/130202823-9a2484d0-c13f-4d86-9685-4bfe04bbf8c2.png" width="90">

[![APHA-CSU](https://circleci.com/gh/APHA-CSU/btb-seq.svg?style=svg)](https://app.circleci.com/pipelines/github/APHA-CSU)

`btb-seq` is the pipeline for APHA's processing of raw *Mycobacterium bovis* Whole Genome Sequencing (WGS) data. The pipeline uses [nextflow](https://www.nextflow.io/docs/latest/getstarted.html) to process batches (1 or more samples) of paired-end `fastq.gz` read files generated on an Illumina sequencer. 

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
- A summary csv file (`FinalOut.csv`) that contains the `Outcome` (see below), WGS Group (Clade) and other high-level metrics for each sample. 
- Consensus `fasta` files
- `.bam` files
- `.vcf` files
- Metagenomics classification of *non-M. Bovis* contaminents

### Run from terminal

To run a batch from the terminal
```
./btb-seq /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/
```

### Run from docker

**Note:** While running from the terminal is the easiest method for developers and data analysts, the pipeline can also be run from docker. This method has the benefit of working across platforms while guaranteeing consistency with automated tests (see below). 

A docker image containing all required dependancies is provided [here](https://hub.docker.com/r/aphacsubot/btb-seq). 

This pull the latest image (if it's not already fetched) from dockerhub and run the container on data
```
sudo docker run --rm -it -v /ABS/PATH/TO/READS/:/reads/ -v /ABS/PATH/TO/RESULTS/:/results/ aphacsubot/btb-seq
```

#### Build docker image from source 

You can also build your own experimental docker image from source
```
docker build ./docker/ -t my-bov-tb
./bov-tb /PATH/TO/READS/ /PATH/TO/OUTPUT/RESULTS/ my-bov-tb
```

## Algorithm

The pipelines processes data in three stages, as shown below. During the preprocessing stage; duplicate reads, low quality bases and adapter sequences are removed from the fastq sample file. Following this, the alignment stage aligns preprocessed reads to a reference genome (*M. bovis* AF2122), performs variant call, masks repeat regions and computes the consensus at each base. The final postprocessing stage assigns an `Outcome` to each sample by analysing data gathered during the preprocessing and alignment stages. The following `Outcome`s are used to signify subsequent lab processing steps:

- **Pass**: The sample contains a known M. Bovis WGS Cluster.
- **Contaminated**: The sample contains contaminants
- **Insufficient Data**: The sample contains insufficient data volumes for sequencing 
- **Check Required**: Further scrutiny of the output is needed as quality thresholds fall below certain criteria, but is likely to contain M.bovis.  

![pipeline](https://user-images.githubusercontent.com/6979169/113730676-ffecef00-96ef-11eb-8670-9fae5e175701.png)


## Validation

This pipeline has been accrediated to ISO17025 standard. It has also been internally validated, tested and approved against a dataset in excess of 10,000 samples that have been sequenced by APHA.  


## Automated Tests

The automated tests provided here ensure the software runs as expected. If you make changes to the algorithm, it is **strongly** reccomended that you run these tests to verify the pipeline is behaving as intended. The tests are also automatically run by `.circleci` on each pull-request. 

### How to run tests

To run a test
```
bash tests/jobs/NAME_OF_TEST.bash
```

### Unit Tests

A number of small tests that assert the functionality of individual components

### Inclusivity Tests

Asserts the `Outcome` and `WGS_CLUSTER` (clade) against samples uploaded by APHA to [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB40340). 

### Limit of Detection (LoD)

The limit of detection test ensures mixtures of M. Avium and M. Bovis at varying proportions give the correct Outcome. This is performed by taking random reads from reference samples of M. Bovis and M. Avium.

| M. Bovis (%) | M. Avium (%) | Outcome |
| ------------- | ------------- | ------------- | 
| 100%   | 0% | Pass | 
| 65%   | 35% | Pass | 
| 60%   | 40% | CheckRequired | 
| 0%   | 100% | Contaminated | 

### Quality Test

The quality test ensures that low quality reads (<20) are not considered for variant calling and genotyping. This is performed by setting uniform quality values to a real-world *M. bovis* sample and asserting output. Low quality bases are removed from the sequence using `Trimmomatic`, which uses a sliding window that deletes reads when the average base quality drops below 20. A table of expected results is shown below.

| Base Quality | Outcome | 
| ------------- | ------------- | 
| 19   | LowQualData | 
| 20   | Pass | 



# Release Process

To release a new version of the software, the `master` branch needs only to be merged into the `prod` branch. To perform this merge, a pull-request from the `master` branch into the `prod` branch needs to be made. Approval of pull-requests to `prod` is made by the CODEOWNER (Richard Ellis). The CODEOWNER is responsible for ensuring the code conforms to the reliability tests. A positive test result is required for approval.

To release a new version of the software:
1. A developer makes a pull-request from the `master` to the `prod` branch. The CODEOWNER is automatically notified by e-mail.
1. The CODEOWNER ensures the automated tests pass on the `master` branch and reviews the code changes. 
1. The CODEOWNER approves the pull-request if they are satisfied, or requests changes.
1. The dev merges the `master` branch into `prod`
1. Following approval, the developer tags the current head of `master` as the next version (see image below). Versions are numbered incrementally with integers, for example `v1`, `v2`, etc. This can be performed by navigating to the github `master` branch and selecting `Create a release`

![image](https://user-images.githubusercontent.com/6979169/163342248-d41c9625-1c79-4463-9425-99522829cd31.png)

![image](https://user-images.githubusercontent.com/6979169/163342279-40bf4673-6af9-4b35-adab-5ea15df601bc.png)

