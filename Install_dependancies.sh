#!/bin/bash
set -e

# This script will fetch and install the dependancies required for the BovTB-nf pipeline.
# The dependancies are generally standard bioinformatics tools
# There are some standard prerequites for a vanilla linux install such as make, make-guile, gcc, zlib-dev, zlib1g-dev,
# libncurses5-dev, libbz2-dev, liblzma-dev, curl, python3, python3-numpy, python3-pip
# e.g. on Ubuntu: sudo apt install make gcc unzip zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev python3 python3-numpy python3-pip
# Followed by pip3 install biopython

# Make directory for dependancy install and cd to that directory before running this script


# FastUniq

wget https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz && tar xzf FastUniq-1.1.tar.gz && rm -f FastUniq-1.1.tar.gz

cd FastUniq/source; make
cd ../..

# Trimmomatic

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && unzip Trimmomatic-0.38.zip && rm -f Trimmomatic-0.38.zip

# bwa

git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..

# samtools and bcftools

wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && tar xjf samtools-1.10.tar.bz2 && rm -f samtools-1.10.tar.bz2
cd samtools-1.10; make 
sudo make install
cd ..

# use this to install latest commit of bcftools (as opposed to the v1.9 release)
#git clone https://github.com/samtools/htslib.git
#git clone https://github.com/samtools/bcftools.git
#cd bcftools; make
#cd ..
 
wget https://github.com/samtools/bcftools/releases/download/1.10.1/bcftools-1.10.1.tar.bz2 && tar xjf bcftools-1.10.1.tar.bz2 && rm -f bcftools-1.10.1.tar.bz2
cd bcftools-1.10.1; make 
sudo make install
cd ..

# bedtools

wget https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz && tar xzf bedtools-2.29.0.tar.gz && rm -f bedtools-2.29.0.tar.gz
cd bedtools2; make
cd ..

# kraken2 and associated database

wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz && tar xzf v2.0.8-beta.tar.gz && rm -f v2.0.8-beta.tar.gz
cd kraken2-2.0.8-beta
./install_kraken2.sh ../Kraken2
cd ..
mkdir Kraken2/db
cd Kraken2/db
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v1_8GB_201904_UPDATE.tgz && tar xvf minikraken2_v1_8GB_201904_UPDATE.tgz && rm -f minikraken2_v1_8GB_201904_UPDATE.tgz
cd ../..

# Add locations to nextflow.config

echo "params.dependPath = "\"$PWD"\"" >> BovTB-nf/nextflow.config
echo "params.kraken2db = "\"$PWD"/Kraken2/db\"" >> BovTB-nf/nextflow.config

echo "Dependancies installed successfully!"

# get some test data

cd ..

mkdir TestData
cd TestData
wget -r -l2 -A ERR84179*.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/ 
mv ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/*/*_1.fastq.gz $PWD
mv ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/*/*_2.fastq.gz $PWD
rm -r ftp.sra.ebi.ac.uk/

cd ..

# anything else??
