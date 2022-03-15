set -e

################## ARGS ##############################

BIOTOOLS_PATH=$PWD/biotools/


################## DEPENDENCIES ######################

apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    make \
    g++ \
    wget
    
################## BIOTOOLS ######################

mkdir -p $BIOTOOLS_PATH
cp ./install*.*sh $BIOTOOLS_PATH
cd $BIOTOOLS_PATH

bash -e install-Kraken2.sh
bash -e install-bracken.sh
