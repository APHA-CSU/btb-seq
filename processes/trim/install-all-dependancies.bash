set -e

################## DEPENDENCIES ######################

sudo apt-get -y update && sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openjdk-11-jdk \
    wget \
    unzip \
    libcurl4-openssl-dev \

################## BIOTOOLS ######################

bash -e install-Trimmomatic.sh
