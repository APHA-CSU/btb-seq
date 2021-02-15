FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Bovine-TB Pipeline Image"
LABEL about.summary="Bioinformatics Pipeline for post-processing of Bovine TB fastq reads"
LABEL about.documentation="https://github.com/APHA-CSU/BovTB-nf"
LABEL about.tags="Genomics, WGS"


################## ARGS #############################

ARG BOVTB_PATH="/BovTB-nf/"
ARG BIOTOOLS_PATH="/biotools/"


################## DEPENDENCIES ######################

# TODO: Remove all except sudo
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo

################## ENTRY ######################

WORKDIR $BOVTB_PATH

# Copy repository
COPY ./ ./
RUN chmod +x ./bin/*

# Install Dependancies
RUN ./docker/build.bash

CMD nextflow run bTB-WGS_process.nf --outdir "/results/" --reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz" -with-report "/results/report.html"
