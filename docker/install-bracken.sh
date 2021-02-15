# bracken
wget https://github.com/jenniferlu717/Bracken/archive/v2.6.0.tar.gz
tar xzf v2.6.0.tar.gz
rm -f v2.6.0.tar.gz
cd Bracken-2.6.0
sh ./install_bracken.sh ../bracken
cd ..
ln -s $PWD/Bracken-2.6.0/bracken /usr/local/bin/bracken