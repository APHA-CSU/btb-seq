set -e

################## ARGS ##############################

BIOTOOLS_PATH=/home/ubuntu/BovTB-nf/biotools/


################## DEPENDENCIES ######################

apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    make \
    g++ \
    wget \
    python3 \
    python3-pip
    
################## BIOTOOLS ######################

mkdir -p $BIOTOOLS_PATH
cp ./install*.*sh $BIOTOOLS_PATH
cd $BIOTOOLS_PATH

# python 
pip3 install biopython numpy pandas gitpython
ln -s /usr/bin/python3 /usr/bin/python

bash -e install-Kraken2.sh
bash -e install-bracken.sh
