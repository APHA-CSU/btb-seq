set -e

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
pip3 install biopython pandas gitpython
sudo ln -s /usr/bin/python3 /usr/bin/python

################## BIOTOOLS ######################


mkdir -p $BIOTOOLS_PATH
cp ./install*.*sh $BIOTOOLS_PATH
cd $BIOTOOLS_PATH

bash -e install-FastUniq.sh
bash -e install-Trimmomatic.sh
bash -e install-bwa.sh
bash -e install-samtools.sh
bash -e install-bcftools.sh
bash -e install-bedtools.sh
bash -e install-Kraken2.sh
bash -e install-bracken.sh
bash -e install-nextflow.sh
bash -e install-sra-toolkit.sh
bash -e install-dwgsim.sh
