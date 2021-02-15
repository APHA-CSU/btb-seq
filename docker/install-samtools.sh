# samtools
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar xjf samtools-1.10.tar.bz2
rm -f samtools-1.10.tar.bz2
cd samtools-1.10
make
sudo make install
cd ..