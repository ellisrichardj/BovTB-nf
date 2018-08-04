#!/bin/bash

	raw_R1=$(zgrep -c "^+" !{forward})	
	uniq_R1=$(grep -c "^+" !{pair_id}_uniq_R1.fastq)
	trim_R1=$(grep -c "^+" !{pair_id}_trim_R1.fastq)
	num_map=$(samtools view -c -F4 !{pair_id}.mapped.sorted.bam)
	avg_depth=$(samtools depth  !{pair_id}.mapped.sorted.bam  |  awk '{sum+=$3} END { print sum/NR}')

	num_raw=$(($raw_R1 * 2))
	num_uniq=$(($uniq_R1 * 2))
	num_trim=$(($trim_R1 * 2))
	pc_aft_dedup=$(echo "scale=2; ($num_uniq*100/$num_raw)" |bc)
	pc_aft_trim=$(echo "scale=2; ($num_trim*100/$num_raw)" |bc)
	pc_mapped=$(echo "scale=2; ($num_map*100/$num_raw)" |bc)

	echo "!{pair_id},"$num_raw","$num_uniq","$pc_aft_dedup","$num_trim","$pc_aft_trim","$num_map","$pc_mapped","$avg_depth"" > !{pair_id}_stats.csv
