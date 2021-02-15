# FastUniq
wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz
tar xzf FastUniq-1.1.tar.gz && rm -f FastUniq-1.1.tar.gz
cd FastUniq/source
make
cd ../..
ln -s $PWD/FastUniq/source/fastuniq /usr/local/bin/fastuniq