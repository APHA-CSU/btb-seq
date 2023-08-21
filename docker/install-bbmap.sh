
wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download -O BBMap_39.01.tar.gz
tar xzf BBMap_39.01.tar.gz && rm BBMap_39.01.tar.gz

sudo ln -s $PWD/bbmap/bbduk.sh /usr/local/bin/bbduk.sh

