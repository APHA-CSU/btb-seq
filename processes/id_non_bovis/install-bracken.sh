# bracken

# TODO: Need to look at alternate solutions running bracken. This works, and is the suggested approach, 
# but there must be a better solution than adding a symlink to the bracken src/ directory

wget https://github.com/jenniferlu717/Bracken/archive/v2.6.0.tar.gz
tar xzf v2.6.0.tar.gz
rm -f v2.6.0.tar.gz
cd Bracken-2.6.0
sh ./install_bracken.sh ../bracken
cd ..
sudo ln -s $PWD/Bracken-2.6.0/bracken /usr/local/bin/bracken
sudo ln -s $PWD/Bracken-2.6.0/src /usr/local/bin/src
  
