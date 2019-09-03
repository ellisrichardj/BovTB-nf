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
*	Version 0.8.2	14/03/19	Define loction of kraken2 database as a nextflow parameter
*	Version 0.8.3	14/03/19	Add option to reduce memory use by kraken2 if required
*	Verison 0.8.4	15/03/19	Correct output location of kraken2 tables
*	Version 0.8.5	25/03/19	Output bam, vcf and consensus to Results directory
*	Version 0.8.6	26/03/19	Add normalization and filtering of indels in vcf before consensus calling
*	Version 0.8.7	29/04/19	Remove intermediary fastq files, rebalance processes and remove redundant process
*	Version 0.8.8	03/05/19	More rebalancing and removing redundant output
*/

params.lowmem = ""
params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"
lowmem = Channel.value("${params.lowmem}")

ref = file(params.ref)
refgbk = file(params.refgbk)
stage1pat = file(params.stage1pat)
stage2pat = file(params.stage2pat)
adapters = file(params.adapters)

pypath = file(params.pypath)
dependpath = file(params.dependPath)
kraken2db = file(params.kraken2db)

/*	Collect pairs of fastq files and infer sample names */
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs } 
	read_pairs.into { read_pairs; raw_reads }


/* remove duplicates from raw data */
process Deduplicate {
	errorStrategy 'ignore'

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

/* trim adapters and low quality bases from fastq data */
process Trim {
	errorStrategy 'ignore'

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

/* map to reference sequence */
process Map2Ref {
	errorStrategy 'ignore'

	publishDir "$params.outdir/Results/bam", mode: 'copy', pattern: '*.sorted.bam'

	maxForks 2

	input:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from trim_read_pairs

	output:
	set pair_id, file("${pair_id}.mapped.sorted.bam") into mapped_bam
	set pair_id, file("${pair_id}.mapped.sorted.bam") into bam4stats

	"""
	$dependpath/bwa/bwa mem -M -t2 $ref  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq |
	 samtools view -@2 -ShuF 2308 - |
	 samtools sort -@2 - -o ${pair_id}.mapped.sorted.bam
	"""
}

/* Variant calling */
process VarCall {
	errorStrategy 'ignore'

	publishDir "$params.outdir/Results/vcf", mode: 'copy', pattern: '*.norm-flt.vcf.gz'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.mapped.sorted.bam") from mapped_bam

	output:
	set pair_id, file("${pair_id}.norm-flt.vcf.gz") into vcf
	set pair_id, file("${pair_id}.norm-flt.vcf.gz") into vcf2
	set pair_id, file("${pair_id}.norm-flt.vcf.gz") into vcf3

	"""
	samtools index ${pair_id}.mapped.sorted.bam
	$dependpath/bcftools/bcftools mpileup -Q 30 -q 60 -Ou -f $ref ${pair_id}.mapped.sorted.bam |
	 $dependpath/bcftools/bcftools call --ploidy 1 -cf GQ - -Ou |
	 $dependpath/bcftools/bcftools norm -f $ref - -Ou |
	 $dependpath/bcftools/bcftools filter --IndelGap 5 –e 'DP<5' -i 'AF>=0.9' - -Oz -o ${pair_id}.norm-flt.vcf.gz
	"""
}

/* Consensus calling */
process VCF2Consensus {
	errorStrategy 'ignore'

	publishDir "$params.outdir/Results/consensus", mode: 'copy', pattern: '*_consensus.fas'
//	publishDir "$params.outdir/Results/bcf", mode: 'copy', pattern: '*.norm-flt.bcf'

	maxForks 2

	input:
	set pair_id, file("${pair_id}.norm-flt.vcf.gz") from vcf2

	output:
	set pair_id, file("${pair_id}_consensus.fas") into consensus //file("${pair_id}.norm-flt.bcf")

	"""
	$dependpath/bcftools/bcftools index ${pair_id}.norm-flt.vcf.gz
#	$dependpath/bcftools/bcftools norm -f $ref ${pair_id}.pileup.vcf.gz -Ob | $dependpath/bcftools/bcftools filter --IndelGap 5 - -Ob -o ${pair_id}.norm-flt.bcf
#	$dependpath/bcftools/bcftools index ${pair_id}.norm-flt.bcf
	$dependpath/bcftools/bcftools consensus -f $ref ${pair_id}.norm-flt.vcf.gz | sed '/^>/ s/.*/>${pair_id}/' - > ${pair_id}_consensus.fas
	"""
}

//	Combine data for generating per sample statistics

raw_reads
	.join(uniq_reads)
	.set { raw_uniq }

trim_reads
	.join(bam4stats)
	.set { trim_bam }

raw_uniq
	.join(trim_bam)
	.set { input4stats }

/* Mapping Statistics*/
process ReadStats{
	errorStrategy 'ignore'

	maxForks 2

	input:
	set pair_id, file("${pair_id}_*_R1_*.fastq.gz"), file("${pair_id}_*_R2_*.fastq.gz"), file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq"), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq"), file("${pair_id}.mapped.sorted.bam") from input4stats

	output:
	set pair_id, file("${pair_id}_stats.csv") into stats
	set pair_id, file('outcome.txt') into Outcome

	shell:
	'''
	raw_R1=$(zgrep -c "^+$" !{pair_id}_*_R1_*.fastq.gz)
	rm !{pair_id}_*_R1_*.fastq.gz
	rm !{pair_id}_*_R2_*.fastq.gz
	uniq_R1=$(grep -c "^+$" !{pair_id}_uniq_R1.fastq)
	rm `readlink !{pair_id}_uniq_R1.fastq`
	rm `readlink !{pair_id}_uniq_R2.fastq`
	trim_R1=$(grep -c "^+$" !{pair_id}_trim_R1.fastq)
	num_map=$(samtools view -c !{pair_id}.mapped.sorted.bam)
	samtools depth -a !{pair_id}.mapped.sorted.bam > depth.txt
	avg_depth=$(awk '{sum+=$3} END { print sum/NR}' depth.txt)
	zero_cov=$(awk '$3<1 {++count} END {print count}' depth.txt)
	sites=$(awk '{++count} END {print count}' depth.txt)
	rm depth.txt

	num_raw=$(($raw_R1*2))
	num_uniq=$(($uniq_R1*2))
	num_trim=$(($trim_R1*2))
	pc_aft_dedup=$(echo "scale=2; ($num_uniq*100/$num_raw)" |bc)
	pc_aft_trim=$(echo "scale=2; ($num_trim*100/$num_raw)" |bc)
	pc_mapped=$(echo "scale=2; ($num_map*100/$num_trim)" |bc)
	genome_cov=$(echo "scale=2; (100-($zero_cov*100/$sites))" |bc)

	mindepth=10
	minpc=60
	minreads=600000
	
	if [ ${avg_depth%%.*} -ge $mindepth ] && [ ${pc_mapped%%.*} -gt $minpc ]; then flag="Pass"
		elif [ ${avg_depth%%.*} -lt $mindepth ] && [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="Comtaminated"
		elif [ ${avg_depth%%.*} -lt $mindepth ] && [ $num_trim -lt $minreads ]; then flag="InsufficientData"
#		elif [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="q_OtherMycobact"
		else flag="CheckRequired"
	fi
 
	echo "Sample,NumRawReads,NumDedupReads,%afterDedup,NumTrimReads,%afterTrim,NumMappedReads,%Mapped,MeanDepth,GenomeCov,Outcome" > !{pair_id}_stats.csv
	echo "!{pair_id},"$num_raw","$num_uniq","$pc_aft_dedup","$num_trim","$pc_aft_trim","$num_map","$pc_mapped","$avg_depth","$genome_cov","$flag"" >> !{pair_id}_stats.csv
	echo "$flag" > outcome.txt
	'''
}

/* SNP filtering and annotation */
process SNPfiltAnnot{

	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.norm-flt.vcf.gz") from vcf3

	output:
	set pair_id, file("${pair_id}.pileup_SN.csv"), file("${pair_id}.pileup_DUO.csv"), file("${pair_id}.pileup_INDEL.csv") into VarTables
	set pair_id, file("${pair_id}.pileup_SN_Annotation.csv") into VarAnnotation

	"""
	$dependpath/bcftools/bcftools view -O v ${pair_id}.norm-flt.vcf.gz | python $pypath/snpsFilter.py - ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp}
	mv _DUO.csv ${pair_id}.pileup_DUO.csv
	mv _INDEL.csv ${pair_id}.pileup_INDEL.csv
	mv _SN.csv ${pair_id}.pileup_SN.csv
	python $pypath/annotateSNPs.py ${pair_id}.pileup_SN.csv $refgbk $ref
	"""
}

//	Combine inputs to assign cluster for each sample

vcf
	.join(stats)
	.set { input4Assign }

/* Assigns cluster by matching patterns of cluster specific SNPs */
process AssignClusterCSS{
	errorStrategy 'ignore'

	maxForks 2

	input:
	set pair_id, file("${pair_id}.norm-flt.vcf.gz"), file("${pair_id}_stats.csv") from input4Assign

	output:
	file("${pair_id}_stage1.csv") into AssignCluster

	"""
	gunzip -c ${pair_id}.norm-flt.vcf.gz > ${pair_id}.pileup.vcf
	python $pypath/Stage1-test.py ${pair_id}_stats.csv ${stage1pat} $ref test 1 ${min_mean_cov} ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} ${min_qual_nonsnp} ${pair_id}.pileup.vcf
	mv _stage1.csv ${pair_id}_stage1.csv
	"""
}

// Collect data to ID any non-M.bovis samples

Outcome
	.join(trim_read_pairs2)
	.set { IDdata }

/* Identify any non-M.bovis samples using kraken */
process IDnonbovis{
	errorStrategy 'ignore'

	publishDir "$params.outdir/Results/NonBovID", mode: 'copy', pattern: '*.tab'

	maxForks 1

	input:
	set pair_id, file('outcome.txt'), file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from IDdata
	val lowmem from lowmem

	output:
	set pair_id, file("${pair_id}_*_kraken2.tab") optional true into IDnonbovis

	"""
	outcome=\$(cat outcome.txt)
	if [ \$outcome != "Pass" ]; then
	$dependpath/Kraken2/kraken2 --threads 2 --quick $lowmem --db $kraken2db --output - --report ${pair_id}_"\$outcome"_kraken2.tab --paired ${pair_id}_trim_R1.fastq  ${pair_id}_trim_R2.fastq 
	else
	echo "ID not required"
	fi
	rm `readlink ${pair_id}_trim_R1.fastq`
	rm `readlink ${pair_id}_trim_R2.fastq`
	"""
}

/* Combine all cluster assignment data into a single results file */
AssignCluster
	.collectFile( name: 'AssignedWGSCluster.csv', sort: true, storeDir: "$PWD/Results", keepHeader: true )



workflow.onComplete {
		log.info "Completed sucessfully:	$workflow.success"		
		log.info "Nextflow Version:	$workflow.nextflow.version"
		log.info "Duration:		$workflow.duration"
		log.info "Output Directory:	$params.outdir/Results"
}
