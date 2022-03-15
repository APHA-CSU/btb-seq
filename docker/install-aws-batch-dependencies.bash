set -e

################## DEPENDENCIES ######################
apt-get update -y
DEBIAN_FRONTEND=noninteractive apt-get install -y awscli 
