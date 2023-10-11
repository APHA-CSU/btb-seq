set -e

################## DEPENDENCIES ######################

sudo apt-get -y update && sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget \
    make \
    liblzma-dev \
    libz-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    gcc \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    g++

# python
sudo ln -s /usr/bin/python3 /usr/bin/python

################## BIOTOOLS ######################

bash -e install-bedtools.sh
bash -e install-bcftools.sh
