# bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz
tar xzf bedtools-2.29.0.tar.gz
rm -f bedtools-2.29.0.tar.gz
cd bedtools2
make
cd ..
ln -s $PWD/bedtools2/bin/bedtools /usr/local/bin/bedtools
