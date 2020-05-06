#!/usr/bin/env nextflow

/*	This is APHA's nextflow pipeline to process Illumina paired-end data from Mycobacterium bovis isolates.  
*	It will first deduplicate the reads using fastuniq, trim them using trimmomatic and then map to the reference genome.
*	Variant positions wrt the reference are determined, togther with data on the number of reads mapping and the depth of 
*	coverage.  Using a panel of predetermined cluster specific SNPs it will also infer cluster membership.
*
*	written by ellisrichardj, based on pipeline developed by Javier Nunez
*
*	Version 0.1.0	31/07/18	Initial version
*	Version 0.2.0	04/08/18	Added SNP filtering and annotation
*	Version 0.3.0	06/08/18	Generate summary of samples in batch
*	Version 0.4.0	16/08/18	Minor improvements
*	Version 0.5.0	14/09/18	Infer genotypes using genotype-specific SNPs (GSS)
*	Version 0.5.1	21/09/18	Fixed bug that prevented GSS from running correctly
*	Version 0.5.2	01/10/18	Mark shorter split hits as secondary in bam file (-M) and change sam flag filter to 3844
*	Version 0.5.3	15/10/18	Increase min mapQ for mpileup to 60 to ensure unique reads only; add $dependpath variable
*	Version 0.6.0	15/11/18	Assigns clusters (newly defined) in place of inferring historical genotype.
*	Version 0.6.1	15/11/18	Fixed bug which caused sample names to be inconsistently transferred between processes
*	Version 0.6.2	24/11/18	Fixed bug to allow cluster assignment to be collected in a single file
*	Version 0.6.3	26/11/18	Removed 'set' in output declaration as it caused nextflow warning
*	Version 0.7.0	26/11/18	Add process to output phylogenetic tree
*	Version 0.7.1	11/12/18	Used join to ensure inputs are properly linked
*	Version 0.7.2	18/12/18	Changed samtools filter to remove unmapped reads and reads that aligned more than once
*	Version 0.7.3	22/01/19	Changed link to adapters file for trimming
*	Version 0.7.4	22/02/19	Changed calculations of mapping stats
*	Version 0.7.5	08/03/19	Further correction of mapping stats and addition of success flags
*	Version 0.8.0	12/03/19	Add kraken to ID samples that fail cluster assignment
*	Version 0.8.1	13/03/19	Update to kraken2
*	Version 0.8.2	14/03/19	Define location of kraken2 database as a nextflow parameter
*	Version 0.8.3	14/03/19	Add option to reduce memory use by kraken2 if required
*	Verison 0.8.4	15/03/19	Correct output location of kraken2 tables
*	Version 0.8.5	25/03/19	Output bam, vcf and consensus to Results directory
*	Version 0.8.6	26/03/19	Add normalization and filtering of indels in vcf before consensus calling
*	Version 0.8.7	29/04/19	Remove intermediary fastq files, rebalance processes and remove redundant process
*	Version 0.8.8	03/05/19	More rebalancing and removing redundant output
*	Version 0.9.0	10/09/19	Filter and mask vcf for consensus calling
*	Version 0.9.1	19/09/19	Remove SNP filtering and annotation process as no longer required
*	Version 0.9.2	20/09/19	Exclude indels from consensus calling step
*	Version 0.9.3   11/10/19	Move ReadStats to standalone shell script
*	Version 0.9.4	17/10/19	Add Data source and datestamp to output files and directories
*	Version 0.9.5	20/12/19	Update to Python3 and lastest samtools/bcftools (1.10)
*	Version 0.9.6	30/03/20	Add Bracken for parsing kraken2 output
*	Version 0.9.7	21/04/20	Determine presence of M.bovis in low quality samples
*	Version 0.9.8	30/04/20	Output snp table in tsv format for each sample
*/

/* Default parameters */
params.lowmem = ""
params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"
lowmem = Channel.value(params.lowmem)

ref = file(params.ref)
refgbk = file(params.refgbk)
rptmask = file(params.rptmask)
stage1pat = file(params.stage1pat)
stage2pat = file(params.stage2pat)
adapters = file(params.adapters)

pypath = file(params.pypath)
dependpath = file(params.dependPath)
kraken2db = file(params.kraken2db)

/*	Collect pairs of fastq files and infer sample names
Define the input raw sequening data files */
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs } 
	read_pairs.into { read_pairs; raw_reads }

// Collect name of data folder and analysis run date
FirstFile = file( params.reads ).first()
	DataPath = FirstFile.getParent()
	GetTopDir = ~/\/\*\//
	TopDir = DataPath - GetTopDir
	params.DataDir = TopDir.last()
	params.today = new Date().format('ddMMMYY')

/* remove duplicates from raw data
This process removes potential duplicate data (sequencing and optical replcaites from the raw data set */
process Deduplicate {
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	set pair_id, file("${pair_id}_*_R1_*.fastq.gz"), file("${pair_id}_*_R2_*.fastq.gz") from read_pairs

	output:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") into dedup_read_pairs
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") into uniq_reads

	"""
	gunzip -c ${pair_id}_*_R1_*.fastq.gz > ${pair_id}_R1.fastq 
	gunzip -c ${pair_id}_*_R2_*.fastq.gz > ${pair_id}_R2.fastq
	echo '${pair_id}_R1.fastq\n${pair_id}_R2.fastq' > fqin.lst
	$dependpath/FastUniq/source/fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
	rm ${pair_id}_R1.fastq
	rm ${pair_id}_R2.fastq
	"""
}	

/* trim adapters and low quality bases from fastq data
Removes the adapters which are added during the lab processing and and any low quality data */
process Trim {
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") from dedup_read_pairs

	output:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_read_pairs
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_read_pairs2
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_reads
	
	"""
	java -jar $dependpath/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 2 -phred33 ${pair_id}_uniq_R1.fastq ${pair_id}_uniq_R2.fastq  ${pair_id}_trim_R1.fastq ${pair_id}_fail1.fastq ${pair_id}_trim_R2.fastq ${pair_id}_fail2.fastq ILLUMINACLIP:$adapters:2:30:10 SLIDINGWINDOW:10:20 MINLEN:36
	rm ${pair_id}_fail1.fastq
	rm ${pair_id}_fail2.fastq
	"""
}

/* map to reference sequence
Aligns the individiual sequence reads to the reference genome */
process Map2Ref {
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/bam", mode: 'copy', pattern: '*.sorted.bam'

	maxForks 2

	input:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from trim_read_pairs

	output:
	set pair_id, file("${pair_id}.mapped.sorted.bam") into mapped_bam
	set pair_id, file("${pair_id}.mapped.sorted.bam") into bam4stats
	set pair_id, file("${pair_id}.mapped.sorted.bam") into bam4mask

	"""
	$dependpath/bwa/bwa mem -M -t2 $ref  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq |
	 $dependpath/samtools-1.10/samtools view -@2 -ShuF 2308 - |
	 $dependpath/samtools-1.10/samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam
	"""
}

/* Variant calling
Determines where the sample differs from the reference genome */
process VarCall {
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/vcf", mode: 'copy', pattern: '*.norm.vcf.gz'

	maxForks 2

	input:
	set pair_id, file("${pair_id}.mapped.sorted.bam") from mapped_bam

	output:
	set pair_id, file("${pair_id}.norm.vcf.gz") into vcf
	set pair_id, file("${pair_id}.norm.vcf.gz") into vcf2

	"""
	$dependpath/samtools-1.10/samtools index ${pair_id}.mapped.sorted.bam
	bcftools mpileup -Q 10 -Ou -f $ref ${pair_id}.mapped.sorted.bam |
	 bcftools call --ploidy 1 -cf GQ - -Ou |
	 bcftools norm -f $ref - -Oz -o ${pair_id}.norm.vcf.gz
	"""
}

/* Masking known repeats regions and sites with zero coverage
Ensure that consensus only includes regions of the genome where there is high confidence */
process Mask {
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	set pair_id, file("${pair_id}.mapped.sorted.bam") from bam4mask

	output:
	set pair_id, file("${pair_id}_RptZeroMask.bed") into maskbed

	"""
	$dependpath/bedtools2/bin/bedtools genomecov -bga -ibam ${pair_id}.mapped.sorted.bam |
	 grep -w "0\$" > ${pair_id}_zerocov.bed
	cat ${pair_id}_zerocov.bed $rptmask | sort -k1,1 -k2,2n |
	 $dependpath/bedtools2/bin/bedtools merge > ${pair_id}_RptZeroMask.bed
	"""
}

// Combine input for consensus calling
// Joins relevant files for the same sample as input to the next process

maskbed
	.join(vcf2)
	.set { vcf_bed }

/* Consensus calling and output snp table for snippy compatability */
process VCF2Consensus {
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/consensus", mode: 'copy', pattern: '*_consensus.fas'
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/snpTables", mode: 'copy', pattern: '*_snps.tab'

	maxForks 2

	input:
	set pair_id, file("${pair_id}_RptZeroMask.bed"), file("${pair_id}.norm.vcf.gz") from vcf_bed

	output:
	set pair_id, file("${pair_id}_consensus.fas") into consensus
	set pair_id, file("${pair_id}_snps.tab") into snpstab

	"""
	bcftools filter --IndelGap 5 -e 'DP<5 && AF<0.8' ${pair_id}.norm.vcf.gz -Ob -o ${pair_id}.norm-flt.bcf
	bcftools index ${pair_id}.norm-flt.bcf
	bcftools consensus -f $ref -e 'TYPE="indel"' -m ${pair_id}_RptZeroMask.bed ${pair_id}.norm-flt.bcf |
	 sed '/^>/ s/.*/>${pair_id}/' > ${pair_id}_consensus.fas
	echo "CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE" > ${pair_id}_snps.tab
	bcftools query -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' ${pair_id}.norm-flt.bcf |
	 awk -F, '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$5":"\$8+\$9" "\$4":"\$6+\$7}' >> ${pair_id}_snps.tab
	"""
}

//	Combine data for generating per sample statistics
// Joins relevant files for the same sample as input to the next process

raw_reads
	.join(uniq_reads)
	.set { raw_uniq }

trim_reads
	.join(bam4stats)
	.set { trim_bam }

raw_uniq
	.join(trim_bam)
	.set { input4stats }

/* Generation of data quality and mapping statistics
Calculate number of raw reads, unique reads, trimmed reads, proportion aligned to reference genome */

process ReadStats{
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 2

	input:
	set pair_id, file("${pair_id}_*_R1_*.fastq.gz"), file("${pair_id}_*_R2_*.fastq.gz"), file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq"), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq"), file("${pair_id}.mapped.sorted.bam") from input4stats

	output:
	set pair_id, file("${pair_id}_stats.csv") into stats
	set pair_id, file('outcome.txt') into Outcome

    """
    ReadStats.sh "$pair_id"
    """
}

//	Combine inputs to assign cluster for each sample
// Joins relevant files for the same sample as input to the next process

vcf
	.join(stats)
	.set { input4Assign }

/* Assigns cluster by matching patterns of cluster specific SNPs
Compares SNPs identified in vcf file to lists in reference table */

process AssignClusterCSS{
	errorStrategy 'finish'
    tag "$pair_id"

	maxForks 1

	input:
	set pair_id, file("${pair_id}.norm.vcf.gz"), file("${pair_id}_stats.csv") from input4Assign

	output:
	file("${pair_id}_stage1.csv") into AssignCluster

	"""
	gunzip -c ${pair_id}.norm.vcf.gz > ${pair_id}.pileup.vcf
	python3 $pypath/Stage1-test.py ${pair_id}_stats.csv ${stage1pat} $ref test ${min_mean_cov} ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} ${min_qual_nonsnp} ${pair_id}.pileup.vcf
	mv _stage1.csv ${pair_id}_stage1.csv
	"""
}

// Collect data to ID any non-M.bovis samples
// Joins relevant files for the same sample as input to the next process

Outcome
	.join(trim_read_pairs2)
	.set { IDdata }

/* Identify any non-M.bovis samples using kraken
Samples with flag != 'Pass' are processed to detemine which microbe(s) are present 
Bracken parses the output which is then sorted to generate a top 20 list of species
Presence / absence of M.bovis is also determined by parsing Bracken output */

process IDnonbovis{
	errorStrategy 'finish'
    tag "$pair_id"

	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/NonBovID", mode: 'copy', pattern: '*.tab'

	maxForks 1

	input:
	set pair_id, file('outcome.txt'), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from IDdata
	val lowmem from lowmem

	output:
	set pair_id, file("${pair_id}_*_brackensort.tab"), file("${pair_id}_*_kraken2.tab")  optional true into IDnonbovis
	file("${pair_id}_bovis.csv") optional true into QueryBovis

	"""
	outcome=\$(cat outcome.txt)
	if [ \$outcome != "Pass" ]; then
	$dependpath/Kraken2/kraken2 --threads 2 --quick $lowmem --db $kraken2db --output - --report ${pair_id}_"\$outcome"_kraken2.tab --paired ${pair_id}_trim_R1.fastq  ${pair_id}_trim_R2.fastq 
	$dependpath/Bracken-2.5.3/bracken -d $kraken2db -r 150 -l S -t 10 -i ${pair_id}_"\$outcome"_kraken2.tab -o ${pair_id}_"\$outcome"_bracken.out
	sed 1d ${pair_id}_"\$outcome"_bracken.out | sort -t \$'\t' -k7,7 -nr - | head -20 > ${pair_id}_"\$outcome"_brackensort.tab
	$dependpath/Bracken-2.5.3/bracken -d $kraken2db -r150 -l S1 -i ${pair_id}_"\$outcome"_kraken2.tab -o ${pair_id}_"\$outcome"-S1_bracken.out
	( sed -u 1q; sort -t \$'\t' -k7,7 -nr ) < ${pair_id}_"\$outcome"-S1_bracken.out > ${pair_id}_"\$outcome"-S1_brackensort.tab
	BovPos=\$(grep 'variant bovis' ${pair_id}_"\$outcome"-S1_brackensort.tab |
	 awk '{print \$1" "\$2" "\$3" "\$4","\$9","(\$10*100)}' || true)
	echo "Sample,ID,TotalReads,Abundance" > ${pair_id}_bovis.csv
	echo "${pair_id},"\$BovPos"" >> ${pair_id}_bovis.csv
	else
	echo "ID not required"
	fi
	rm `readlink ${pair_id}_trim_R1.fastq`
	rm `readlink ${pair_id}_trim_R2.fastq`
	"""
}

/* Combine all cluster assignment data into a single results file */

AssignCluster
	.collectFile( name: "${params.DataDir}_AssignedWGSCluster_${params.today}.csv", sort: true, storeDir: "$PWD/Results_${params.DataDir}_${params.today}", keepHeader: true )

QueryBovis
	.collectFile( name: "${params.DataDir}_BovPos_${params.today}.csv", sort: true, storeDir: "$PWD/Results_${params.DataDir}_${params.today}", keepHeader: true )

workflow.onComplete {
		log.info "Completed sucessfully:	$workflow.success"		
		log.info "Nextflow Version:	$workflow.nextflow.version"
		log.info "Duration:		$workflow.duration"
		log.info "Output Directory:	$params.outdir/Results_${params.DataDir}_${params.today}"
}
