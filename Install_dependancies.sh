#!/bin/bash
set -e

# This script will fetch and install the dependancies required for the BovTB-nf process.
# The dependancies are generally standard bioinformatics tools
# There are some standard prerequites for a vanilla linux install such as make, make-guile, gcc, zlib-dev, zlib1g-dev,
# libncurses5-dev, libbz2-dev, liblzma-dev, python (not python3), python-numpy, python-pip
# e.g. on Ubuntu: sudo apt install make make-guile gcc zlib-dev zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev python python-numpy, python-pip
# Followed by pip install biopython


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

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xjf samtools-1.9.tar.bz2 && rm -f samtools-1.9.tar.bz2
cd samtools-1.9; make 
sudo make install
cd ..
 
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && tar xjf bcftools-1.9.tar.bz2 && rm -f bcftools-1.9.tar.bz2
cd bcftools-1.9; make 
sudo make install
cd ..

# vcfutils.pl

wget https://github.com/lh3/samtools/blob/master/bcftools/vcfutils.pl

# kraken2 and associated database

wget http://github.com/DerrickWood/kraken2/archive/v2.0.7-beta.tar.gz && tar xzf v2.0.7-beta.tar.gz && rm -f v2.0.7-beta.tar.gz
cd kraken2-2.0.7-beta
./install_kraken2.sh ../Kraken2
cd ..
mkdir Kraken2/db
cd Kraken2/db
wget https://ccb.jhu.edu/software/kraken2/dl/minikraken2_v1_8GB.tgz && tar xvf minikraken2_v1_8GB.tgz
cd ../..
export KRAKEN2_DEFAULT_DB="$PWD/Kraken2/db/:"

# get some test data

#mkdir Data
#cd Data
#wget -r -l2 -A ERR84179*.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/ 
#mv ftp.sra.ebi.ac.uk/vol1/fastq/ERR841//*_1.fastq.gz $PWD
#mv ftp.sra.ebi.ac.uk/vol1/fastq/ERR841//*_2.fastq.gz $PWD
#rm -r ftp.sra.ebi.ac.uk/

# anything else??



# Add locations to nextflow.config

echo "params.dependPath = "\"$PWD"\"" >> BovTB-nf/nextflow.config

echo "Dependancies installed sucessfully!"

