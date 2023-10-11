# bwa
git clone https://github.com/lh3/bwa.git
cd bwa
git checkout v0.7.17
make
cd ..
sudo ln -s $PWD/bwa/bwa /usr/local/bin/bwa
