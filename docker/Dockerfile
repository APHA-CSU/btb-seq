FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Bovine-TB Pipeline Image"
LABEL about.summary="Bioinformatics Pipeline for post-processing of Bovine TB fastq reads"
LABEL about.documentation="https://github.com/APHA-CSU/BovTB-nf"
LABEL about.tags="Genomics, WGS"


################## ARGS #############################

ARG BOVTB_PATH="/home/ubuntu/BovTB-nf/"

################## DEPENDENCIES ######################

# Copy repository
WORKDIR $BOVTB_PATH
COPY ./ ./
RUN chmod +x ./bin/*

# Sudo
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive

# Biotools
RUN cd ./docker && bash ./install-all-dependancies.bash

## AWS Batch
RUN mkdir ./nxf-scratch-dir
RUN cd ./docker && bash ./install-aws-batch-dependencies.bash
