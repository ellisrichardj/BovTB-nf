#**BovTB-nf**

------------

This is the updated pipeline for APHA's processing of *Mycobacterium bovis* WGS data. BovTB-nf is designed to process a batch (1 or more samples) of paired-end fastq files generated on an Illumina sequencer. It will first remove duplicate reads from the dataset (FastUniq) and then trim the unique reads based on base-call quality and the presence of adapters (Trimmomatic). Reads are then mapped to the *M. bovis* AF2122 reference genome and variants called (bwa/samtools/bcftools).

It has been built to run using nextflow, using standard bioinformatic tools for the most part. The external dependancies are:
-	FastUniq
-	Trimmomatic
-	bwa
-	samtools and bcftools
-	vcfutils.pl

##Installation

Of course Nextflow itself is a prerequisite and should be installed as described in the [Nextflow Documentation](https://www.nextflow.io/docs/latest/getstarted.html)

If you have the dependancies installed the pipeline can run by simply typing: 

	nextflow run ellisrichardj/BovTB-nf

Alternatively, clone the repository:

	git clone https://github.com/ellisrichardj/BovTB-nf.git

If required, there is simple script for installing the dependancies (helpfully called Install_dependancies.sh), which will also update the nextflow config file with their locations.

-------------

##Examples

In its simplest form just run the Nextflow process from the directory containing the fastq files:

	cd /path/to/Data
	nextflow run BovTB-nf




