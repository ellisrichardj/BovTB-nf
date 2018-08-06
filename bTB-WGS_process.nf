#!/usr/bin/env nextflow

/*	This is nextflow pipeline to process Illumina paired-end data from Mycobacterium bovis isolates  
*	It will first deduplicate the reads using fastuniq, trim them using trimmomatic and then map to the reference genome.
*	Variant positions wrt the reference are determined, togther with data on the number of reads mapping and the depth of 
*	coverage.  
*
*	written by ellisrichardj, based on pipeline developed by Javier Nunez
*
*	Version 0.1.0	31/07/18	Initial version
*	Version 0.2.0	04/08/18	Added SNP filtering and annotation
*	Version 0.3.0	06/08/18	Generate summary of samples in batch
*/


params.reads = "$PWD/*_{S*_R1,S*_R2}*.fastq.gz"
params.outdir = "$PWD"

ref = file(params.ref)
refgbk = file(params.refgbk)
stage1pat = file(params.stage1pat)
stage2pat = file(params.stage2pat)

pypath = file(params.pypath)


/*	Collect pairs of fastq files and infer sample names */
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs} 
	read_pairs.into { read_pairs; raw_reads }


/* remove duplicates from raw data */
process Deduplicate {
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file(forward), file(reverse) from read_pairs

	output:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") into dedup_read_pairs
	file("${pair_id}_uniq_R1.fastq") into uniq_reads

	"""
	gunzip -c ${forward} > ${pair_id}_R1.fastq 
	gunzip -c ${reverse} > ${pair_id}_R2.fastq
	echo '${pair_id}_R1.fastq\n${pair_id}_R2.fastq' > fqin.lst
	~/FastUniq/source/fastuniq -i fqin.lst -o ${pair_id}_uniq_R1.fastq -p ${pair_id}_uniq_R2.fastq
	"""
}	

/* trim adapters and low quality bases from fastq data */
process Trim {
	errorStrategy 'ignore'

	maxForks 1

	input:
	set pair_id, file("${pair_id}_uniq_R1.fastq"), file("${pair_id}_uniq_R2.fastq") from dedup_read_pairs

	output:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") into trim_read_pairs
	file("${pair_id}_trim_R1.fastq") into trim_reads
	
	"""
	java -jar ~/MyScripts/trimmomatic-0.30.jar PE -threads 4 -phred33 ${pair_id}_uniq_R1.fastq ${pair_id}_uniq_R2.fastq  ${pair_id}_trim_R1.fastq ${pair_id}_fail1.fastq ${pair_id}_trim_R2.fastq ${pair_id}_fail2.fastq ILLUMINACLIP:/home/richard/ReferenceSequences/adapter.fasta:2:30:10 SLIDINGWINDOW:10:20 MINLEN:36
	rm ${pair_id}_fail1.fastq
	rm ${pair_id}_fail2.fastq
	"""
}

/* map to reference sequence */
process Map2Ref {
	errorStrategy 'ignore'

	maxForks 1

	input:
	set pair_id, file("${pair_id}_trim_R1.fastq"), file("${pair_id}_trim_R2.fastq") from trim_read_pairs

	output:
	set pair_id, file("${pair_id}.mapped.sorted.bam") into mapped_bam
	file("${pair_id}.mapped.sorted.bam") into bam4stats

	"""
	bwa mem -T10 -t4 $ref  ${pair_id}_trim_R1.fastq ${pair_id}_trim_R2.fastq |
	 samtools view -@4 -ShuF 4 - |
	 samtools sort -@4 - -o ${pair_id}.mapped.sorted.bam
	"""
}

/* Variant calling */
process VarCall {
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.mapped.sorted.bam") from mapped_bam

	output:
	set pair_id, file("${pair_id}.pileup.vcf.gz") into vcf, vcf2

	"""
	samtools index ${pair_id}.mapped.sorted.bam
	samtools mpileup -q 10 -uvf $ref ${pair_id}.mapped.sorted.bam |
	 bcftools call --ploidy Y -cf GQ - -Oz -o ${pair_id}.pileup.vcf.gz
	"""
}

/* Consensus calling */
process VCF2Consensus {
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.pileup.vcf.gz") from vcf

	output:
	set pair_id, file("${pair_id}_consensus.fas") into consensus

	"""
	bcftools index ${pair_id}.pileup.vcf.gz
	bcftools consensus -f $ref -o ${pair_id}_consensus.fas ${pair_id}.pileup.vcf.gz
	"""
}


/* Mapping Statistics*/
process ReadStats{
	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file(forward) from raw_reads
	file("${pair_id}_uniq_R1.fastq") from uniq_reads
	file("${pair_id}_trim_R1.fastq") from trim_reads
	file("${pair_id}.mapped.sorted.bam") from bam4stats

	output:
	file("${pair_id}_stats.csv") into stats

	shell:
	'''
	raw_R1=$(zgrep -c "^+" !{forward})	
	uniq_R1=$(grep -c "^+" !{pair_id}_uniq_R1.fastq)
	trim_R1=$(grep -c "^+" !{pair_id}_trim_R1.fastq)
	num_map=$(samtools view -c -F4 !{pair_id}.mapped.sorted.bam)
	avg_depth=$(samtools depth  !{pair_id}.mapped.sorted.bam  |  awk '{sum+=$3} END { print sum/NR}')

	num_raw=$(($raw_R1*2))
	num_uniq=$(($uniq_R1*2))
	num_trim=$(($trim_R1*2))
	pc_aft_dedup=$(echo "scale=2; ($num_uniq*100/$num_raw)" |bc)
	pc_aft_trim=$(echo "scale=2; ($num_trim*100/$num_raw)" |bc)
	pc_mapped=$(echo "scale=2; ($num_map*100/$num_raw)" |bc)

	echo "!{pair_id},"$num_raw","$num_uniq","$pc_aft_dedup","$num_trim","$pc_aft_trim","$num_map","$pc_mapped","$avg_depth"" > !{pair_id}_stats.csv
	'''
}

/* SNP filtering and annotation */
process SNPfiltAnnot{

	errorStrategy 'ignore'

	maxForks 4

	input:
	set pair_id, file("${pair_id}.pileup.vcf.gz") from vcf2

	output:
	set pair_id, file("${pair_id}.pileup_SN.csv"), file("${pair_id}.pileup_DUO.csv"), file("${pair_id}.pileup_INDEL.csv"), file("${paid_id}.pileup_Annotation.csv") into VarTables
	set pair_id, file("${pair_id}.pileup_SN_Annotation.csv") into VarAnnotation

	"""
	gunzip -c ${pair_id}.pileup.vcf.gz > ${pair_id}.pileup.vcf
	python $pypath/snpsFilter.py ${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} ${pair_id}.pileup.vcf
	python $pypath/annotateSNPs.py ${pair_id}.pileup_SN.csv $refgbk $ref
	rm ${pair_id}.pileup.vcf
	"""
}

/* Genotyping - inference of spoligo and VNTR types
process Genotyping{

	input:

	output:


	"""
**variables need defining**

	python Stage1_TBRun_2.py 
**$results_path** 
$stage1pat $ref 
**$run_name** 
**$str(ncores_stage2)** 
**$str(th_mean_coverage_mapping)** 
${min_cov_snp} ${alt_prop_snp} ${min_qual_snp} 
**$str(th_quality_nonSNP)**

	python Stage2_TBRun_2.py 
**$data_path** 
**$results_path** 
$stage2pat 
**$run_name** 
**$str(ncores_stage2)** 
**$str(th_mean_coverage_mapping)** 
**$str(stage1_done).lower()** 
**$str(th_line_contamination)** 
**$str(stats_test)**

}
*/

/* Run summary */
process CombineRun{

	errorStrategy 'ignore'

	input:
	file("${pair_id}_stats.csv") from stats

	output:
	file("RunStats.csv") into SummaryStats

	"""
	cat *_stats.csv > RunStats.csv
	"""

}

workflow.onComplete {
		log.info "Completed sucessfully:	$workflow.success"		
		log.info "Nextflow Version:	$workflow.nextflow.version"
		log.info "Duration:		$workflow.duration"
		log.info "Output Directory:	$params.outdir"
}
