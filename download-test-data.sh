# This file downloads some test data for the Bovine TB pipeline

cd ..

mkdir TestData
cd TestData
wget -r -l2 -A ERR84179*.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/ 
mv ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/*/*_1.fastq.gz $PWD
mv ftp.sra.ebi.ac.uk/vol1/fastq/ERR841/*/*_2.fastq.gz $PWD
rm -r ftp.sra.ebi.ac.uk/

cd ..
