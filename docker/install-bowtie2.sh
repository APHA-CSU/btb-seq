#!/bin/bash

wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.1.tar.gz
tar xzf v2.5.1.tar.gz && rm -f v2.5.1.tar.gz 
cd bowtie2-2.5.1/
make
sudo ln -s ~/biotools/bowtie2-2.5.1/bowtie2 /usr/local/bin/bowtie2
