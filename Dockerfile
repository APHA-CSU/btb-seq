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

RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openjdk-11-jdk \
    sudo \
    wget \
    make \
    git \
    curl \
    liblzma-dev \
    libz-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    gcc \
    unzip \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-numpy \
    python3-pip \
    && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/*

# python 
RUN pip3 install biopython pandas && \
    ln -s /usr/bin/python3 /usr/bin/python


################## BIOTOOLS ######################

WORKDIR $BIOTOOLS_PATH

# FastUniq
RUN wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz && \
    tar xzf FastUniq-1.1.tar.gz && rm -f FastUniq-1.1.tar.gz && \
    cd FastUniq/source && \
    make && \
    cd ../..

# Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && \
    unzip Trimmomatic-0.38.zip && \
    rm -f Trimmomatic-0.38.zip

# bwa
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    cd ..

# samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar xjf samtools-1.10.tar.bz2 && \
    rm -f samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    make && \
    sudo make install && \
    cd ..

# use this to install latest commit of bcftools (as opposed to the v1.9 release)
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
    tar xjf bcftools-1.10.2.tar.bz2 && \
    rm -f bcftools-1.10.2.tar.bz2 && \
    cd bcftools-1.10.2 && \
    make && \
    sudo make install && \
    cd ..

# bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz && \
    tar xzf bedtools-2.29.0.tar.gz && \
    rm -f bedtools-2.29.0.tar.gz && \
    cd bedtools2 && \
    make && \
    cd ..

# kraken2 and associated database
RUN wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz && \
    tar xzf v2.0.8-beta.tar.gz && \
    rm -f v2.0.8-beta.tar.gz && \
    cd kraken2-2.0.8-beta && \
    ./install_kraken2.sh ../Kraken2 && \
    cd ..

RUN mkdir Kraken2/db && \
    cd Kraken2/db && \
    wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz && \
    tar xvf minikraken2_v1_8GB_201904.tgz && \
    rm -f minikraken2_v1_8GB_201904.tgz && \
    cd ../..

# bracken
RUN wget https://github.com/jenniferlu717/Bracken/archive/v2.5.3.tar.gz && \
    tar xzf v2.5.3.tar.gz && \
    rm -f v2.5.3.tar.gz && \
    cd Bracken-2.5.3 && \
    sh ./install_bracken.sh ../bracken && \
    cd ..

# Install nextflow.
# ENV line required to set jvm memory and cpu limits to docker. 
# see: https://github.com/nextflow-io/nextflow/blob/v20.07.1/docker/Dockerfile
COPY install_nextflow-20.7.1.bash ./install_nextflow.bash
RUN cat ./install_nextflow.bash | bash
RUN ln -s $PWD/nextflow /usr/local/bin/nextflow

################## ENTRY ######################

WORKDIR $BOVTB_PATH

# Copy repository
COPY ./ ./

# Add locations to nextflow.config
RUN echo "params.dependPath = \"$BIOTOOLS_PATH\"" >> ./nextflow.config
RUN echo "params.kraken2db = \"$BIOTOOLS_PATH/Kraken2/db/minikraken2_v1_8GB/\"" >> ./nextflow.config

CMD /bin/bash
