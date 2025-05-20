FROM ubuntu:24.04

################## METADATA ##########################

LABEL base.image=ubuntu:24.04
LABEL software="Bovine-TB Pipeline Image"
LABEL about.summary="Bioinformatics Pipeline for post-processing of Bovine TB fastq reads"
LABEL about.documentation="https://github.com/APHA-CSU/btb-seq"
LABEL about.tags="Genomics, WGS"

################## DEPENDENCIES ######################

RUN apt-get -y update && apt-get install --yes --no-install-recommends \
    wget \
    curl \
    libcurl4-openssl-dev \
    git \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    build-essential \
    default-jre \
    nano \
    bc \
    unzip

RUN apt-get -y update && apt-get install --yes --no-install-recommends \
    python3 \
    python3-pip \
    python3-pandas \
    python3-numpy \
    python3-venv \
    python3-biopython \
    python3-git \
    python-is-python3

RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 --no-check-certificate && \
    tar xjf samtools-1.21.tar.bz2 && \
    rm -f samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    make && \
    make install && \
    cd ..

RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 --no-check-certificate && \
    tar xjf bcftools-1.21.tar.bz2 && \
    rm -f bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    make && \
    make install && \
    cd ..

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz --no-check-certificate && \
    tar xzf bedtools-2.31.1.tar.gz && \
    rm -f bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    ln -s $PWD/bin/bedtools /usr/local/bin/bedtools && \
    cd ..

RUN wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz --no-check-certificate && \
    tar xzf v2.0.8-beta.tar.gz && \
    rm -f v2.0.8-beta.tar.gz && \
    cd kraken2-2.0.8-beta && \
    ./install_kraken2.sh ../Kraken2 && \
    cd ..  && \
    ln -s $PWD/Kraken2/kraken2 /usr/local/bin/kraken2

RUN mkdir -p /opt/Kraken2/db  && \
    wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v1_8GB_201904.tgz -nv --no-check-certificate && \
    tar xvf minikraken2_v1_8GB_201904.tgz -C /opt/Kraken2/db/  && \
    rm -f minikraken2_v1_8GB_201904.tgz

RUN wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v3.1.tar.gz --no-check-certificate && \
    tar xvf v3.1.tar.gz && \
    rm -f v3.1.tar.gz && \
    cd Bracken-3.1/ && \
    sh ./install_bracken.sh && \
    cd .. && \
    ln -s $PWD/Bracken-3.1/bracken /usr/local/bin/bracken

RUN wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.19.tar.gz --no-check-certificate && \
    tar xvf v0.7.19.tar.gz && \
    rm -f v0.7.19.tar.gz && \
    cd bwa-0.7.19 \
    make && \
    cd .. && \
    ln -s $PWD/bwa-0.7.19/bwa /usr/local/bin/bwa

RUN wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz --no-check-certificate && \
    tar xzf FastUniq-1.1.tar.gz && rm -f FastUniq-1.1.tar.gz && \
    cd FastUniq/source && \
    make && \
    cd ../.. && \
    ln -s $PWD/FastUniq/source/fastuniq /usr/local/bin/fastuniq

RUN wget https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-0.39.zip --no-check-certificate && \
    unzip Trimmomatic-0.39.zip && rm Trimmomatic-0.39.zip && \
    ln -s $PWD/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar
