BIOTOOLS_PATH="~/biotools/"


################## DEPENDENCIES ######################

sudo apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openjdk-11-jdk \
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
    vim \
    nano \
    bc

# python 
pip3 install biopython pandas
ln -s /usr/bin/python3 /usr/bin/python


################## BIOTOOLS ######################

mkdir -p $BIOTOOLS_PATH
cd $BIOTOOLS_PATH

sh install-FastUniq.sh
sh install-Trimmomatic.sh
sh install-bwa.sh
sh install-samtools.sh
sh install-bcftools.sh
sh install-bedtools.sh
sh install-Kraken2.sh
sh install-bracken.sh
sh install-nextflow.sh
sh install-sra-toolkit.sh

#### TODO: is this good?
chmod +x ./bin/*
