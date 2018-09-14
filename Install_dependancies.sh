#!/bin/bash
set -e

# This script will fetch and install the dependancies required for the BovTB-nf process.
# The dependancies are generally standard bioinformatics tools


# FastUniq

wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz && tar xzf FastUniq-1.1.tar.gz && rm -f FastUniq-1.1.tar.gz

cd FastUniq/source; make
cd ../..

# Trimmomatic

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && unzip Trimmomatic-0.38.zip && rm -f Trimmomatic-0.38.zip

# bwa

git clone https://github.com/lh3/bwa.git
cd bwa; make

# samtools and bcftools

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xjf samtools-1.9.tar.bz2 && rm -f samtools-1.9.tar.bz2
cd samtools-1.9; make 
make install
cd ..
 
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar xjf bcftools-1.9.tar.bz2 && rm -f bcftools-1.9.tar.bz2
cd bcftools-1.9; make 
make install
cd ..

# vcfutils.pl

wget https://github.com/lh3/samtools/blob/master/bcftools/vcfutils.pl

# anything else??



# Add locations to nextflow.config

echo "params.dependPath = "$PWD"" > nextflow.config
