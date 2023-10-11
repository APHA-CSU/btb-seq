set -e

################## DEPENDENCIES ######################

sudo apt-get -y update && sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    git \
    python3 \
    python3-pip \

# python 
pip3 install biopython numpy pandas gitpython
sudo ln -s /usr/bin/python3 /usr/bin/python
