################## ARGS ##############################

BIOTOOLS_PATH=~/biotools/


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
sudo ln -s /usr/bin/python3 /usr/bin/python

################## BIOTOOLS ######################

echo BIOTOOLS $BIOTOOLS_PATH

mkdir -p $BIOTOOLS_PATH

echo mkdir

cp ./install*.sh $BIOTOOLS_PATH

echo cp

cd $BIOTOOLS_PATH

echo cd


sh -e install-FastUniq.sh
sh -e install-Trimmomatic.sh
sh -e install-bwa.sh
sh -e install-samtools.sh
sh -e install-bcftools.sh
sh -e install-bedtools.sh
sh -e install-Kraken2.sh
sh -e install-bracken.sh
sh -e install-nextflow.sh
sh -e install-sra-toolkit.sh
