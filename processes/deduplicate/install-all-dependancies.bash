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
    unzip \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \

################## BIOTOOLS ######################

bash -e install-FastUniq.sh
