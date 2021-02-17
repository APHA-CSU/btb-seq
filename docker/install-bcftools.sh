# use this to install latest commit of bcftools (as opposed to the v1.9 release)
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xjf bcftools-1.10.2.tar.bz2
rm -f bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2
make
sudo make install
cd ..