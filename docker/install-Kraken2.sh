# kraken2 and associated database
wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
tar xzf v2.0.8-beta.tar.gz
rm -f v2.0.8-beta.tar.gz
cd kraken2-2.0.8-beta
./install_kraken2.sh ../Kraken2
cd ..
ln -s $PWD/Kraken2/kraken2 /usr/local/bin/kraken2

mkdir Kraken2/db

wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v1_8GB_201904.tgz
tar xvf minikraken2_v1_8GB_201904.tgz -C Kraken2/db/
rm -f minikraken2_v1_8GB_201904.tgz