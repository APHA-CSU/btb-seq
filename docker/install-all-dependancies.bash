set -e

################## ARGS ##############################

BIOTOOLS_PATH=biotools/


################## DEPENDENCIES ######################

apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    make \
    wget
    
################## BIOTOOLS ######################

mkdir -p $BIOTOOLS_PATH
cp ./install*.*sh $BIOTOOLS_PATH
cd $BIOTOOLS_PATH

bash -e install-Kraken2.sh
bash -e install-bracken.sh
