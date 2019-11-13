#!/usr/bin/python
# python convert_Mobster_predictions_into_Refined_predictions.py -i input_Mobster -b input_bam_reads -l read_length_of_bam_file_reads -r input_reference_genome_fasta -o output_Mobster [-s restart_input_Mobster] [-q quality_filter] [-d discard_Mobster_MEIs] [-c num_cores] [-k keep_chr_prefix]
# python convert_Mobster_predictions_into_Refined_predictions.py \
# 	-i MGRBp1_sample357_chr11_fourMEI_includeMEregions_predictions.txt \
# 	-b AAARX_chr11.bam \
# 	-l 150 \
# 	-r ~/projects_associated_with_mobile_elements_and_telomeres_and_mitochondrial_DNA/reference_genomes/hs37d5x/hs37d5x.fa \
# 	-e ~/projects_associated_with_mobile_elements_and_telomeres_and_mitochondrial_DNA/mobile_elements_investigations/a_Mobster_environment/Mobster_environment/Mobster/repmask/hg19_alul1svaerv.txt \
# 	-o MGRBp1_sample357_chr11_fourMEI_includeMEregions_predictions_refined_afterChange2.txt \
# 	-q 16 \
# 	-d MGRBp1_sample357_chr11_fourMEI_discarded_predictions.txt \
#	-c 16

# In this program:
#     * MEI is mobile element insertion
#     * border5 is the position in the genome that contains the ref.seq., and 1 nucleotide to the left of it is the 3-prime end of the MEI
#     * border3 is the position in the genome that contains the ref.seq., and 1 nucleotide to the right of it is the 5-prime end of the MEI

# Input file "input_Mobster" contains multiple mobile element insertions (MEI), 1 per line, for one sample.
# For the MEI that Mobster identified by split reads, the border5 and border3 (3-prime and 5-prime ends of the mobile element respectively) are often exact.
# Sometimes for MEI identified by split reads, border5 and border3 are not exact. 
# For MEI identified by discordant read pairs (one read mapped to one place and the read-pair mapped far away from that)
# the border5 and border3 are not exact and can't be exact because the exact insertion points have not been read and seen.
# Thus, for a given Mobster-identified MEI, we don't know if it has identified exact insertion points or not.

# For each MEI, this program obtains the bam reads in that small region and the reference sequence for that region.
# The program then tries to identify the exact insertion points.
# If it can identify the exact insertion points, then those points and the MEI sequences are appended to the end of the Mobster record.
# The MEI sequences are those ends of sequences that don't match the reference genome and are truncated because they no longer match the reference genome.

# The input_reference_genome_fasta needs to have been indexed with samtools faidx.

# If the restart_input_Mobster_file input parameter is present, then this program is restarted from that point.
# The restart_input_Mobster_file is copied to output and the input file records are skipped until after that point, from then on processing continues.

# If the quality filter input parameter is present, then MEI calls having quality less than that are discarded.
# The quality is calculated using Mobster fields: 
#	cluster5_hits, cluster3_hits, split5_hits, split3_hits, 
#	polyA5_hits, polyT5_hits, polyA3_hits, polyT3_hits, original discordant unique.
# PHRED-like scale: ie. 30 = 1 in 1000 chance of not being a ME, 20 = 1 in 100, 10 = 1 in 10 chance.
#	= log( (2*split5 + 2*split3 + discord5 + discord3 + polyA5+polyT5+polyA3+polyT3 + orig_discord_unique)*10 )

# If the discarded_MEI_calls_file input parameter is present, then discarded MEI calls are written to this file.
# MEI calls are discarded when this program finds no support for their existance, 
# when their quality score is lower than the filter quality score,
# or when they are mapped to a chromosome or contig to be ignored.

# If the input parameter -k for keep_chr_prefix, then don't remove the chr in front of chromosomes.
# If the value is SOME, then don't remove the chr for 1-22,X,Y,MT, but do remove it for everything else, such as chrGL000193 becomes GL000193

# This program obtains the bams reads from identify_split_read_clusters.py.

# For are examples of bams reads returned from identify_split_read_clusters.py for 4 different regions of a given sample
# showing that identify_split_read_clusters.py sometimes shows multiple sequences that are almost the same
# and need to be considered the same by this program.
# Example has one mobile elemement insertions, giving 2 insertion points (for the 3-prime and 5-prime ends of the mobile element).
#
# CTTCTCTCTCCTCCTCCTCCGCCAGGTACACAGCCTTCTNATTTGGCATCCCAGACACCACAGCAACAGATGGGGCATCGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATC                                  
# CTTCTCTCTCCTCCTCCTCCGCCAGGTACACAGCCTTCTTATTTGGCATCCCAGACACCACAGCAACAGATGGGGCATCTTGTCTCGCATGGCTGTGGGATATCTGACTGGCCAAGGCCAGGACTCTGTTGAGGGGGTAGGGATGCCAGCGTTGTCCACGGCACTACAGATTCT
# CTTCTCTCTCCTCCTCCTCCGCCAGGTACACAGCCTTCTNATTTGGCATCCCAGACACCACAGCAACAGATGGGGCATCGCCGGGCGCGGTGCTGTGGGATATCTGACTGGCCAAGGCAAGGACCCTGTTGAGGGGGTAGGGATGCCAGCTTTGTCCACGGCACTACAGATTCT
# CTTCTCCCTCCTCCTCCTCCGCCAGGTACACAGCCTTCTTATTTGGCATCCCAGACACCACAGCAACAGATGGGGCATCTTGTCTCGCATGGCTGTGGGATATCTGACTGGCCAAGGCAAGGACTCTGTTGAGGGGGTAGGGATGCCAGCTTTGTCCACGGCACTACAGATTCT
#          GACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAACAGATGGGGCATCTTGTCTCGCATGGCTGTGGGATATCTGACTGGCCAAGGCAAGGACTCTGTTGAGGGGGTAGGGATGCCAGCCTTGTCCACGGCACTACAGATTCT
#
# CACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTGTTCCTCAAATGAGTAAAAGGCACTTCTTTACTGCTGGAATGAGGCATGTAGGTGAAAGTAAGTTGAAGTAGTCGAATAACATCTATACAGTGAGTCCTGCAAGACTT
# CACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTATTCGTAAAATGAGTAAAAGGCACTTCTGTACTGCTGGAATGAGGCATGTAGTTGAAAGTAAGTTGAAGTAGTCGAATAACATCGATCCAGTGAGTCCTGCAAGACTT
# CACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTATTCCTGAAAAGAGTAAAAGGCACTTCTGTACTGCTGGAATGAGGCATGGAGTGGAAAGGAAGTTGAAGTAGTCTGATAACAACTATCCAGGG               
# CACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTGTTCCTCAAATGAGTAAAAGGCACTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGGAGCGGAGGTCCTCGCT                                     
# CACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTGTTCCTCAAATGAGTAAAAGGCACTTCTTTACTGCTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGAATAACATCTATCCCGTGAGTCCTGCAAGACTT
#      ATGATCCACCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGAAATGAGTAAAAGGCACTTCTTTACTGCTGGAATGAGGCATGTAGGTGAAAGTAAGTTGAAGTAGTCGAATAACATCTATACAGTGAGTCCTGCAAGACTT
# 
# GTGTGTTTTTGTTGTTTTGACCCTTCTCACGTTGTTTTTTAGTATGCCCCATTTGTCCACAATTTTAACAGGACCCTGAGAACTTTTGTGTGTTATTTCCAACTTTTAAACCTGCCATTGCCGGGCCTAGGAGGGTAGCTTTGTATACATGTCCTCCAATG
# GTGTGTTTTTGTTGTTTTGACCCTTCTCACGTTGTTTTTTAGTATGCCCCATTTGTCCACAATTTTAACAGGACCCTGAGAACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGCGGGGGCTC                                   
# GTGGGTTTTTGTTGTTTTGACCCTTCTCACGTTGTTTTTTAGTATGCCCCATTTGTCCACAATTGTAACAGGACCCTGAGAACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTAGGGGGGGTTCCCGCTTGGTCCCCTGCTTG                 
# GTGTGTTTTTGTTGTTTTGACCCTTCTCACGTTGTTTTTTAGTATGCCCCATTTGTCCACAATTTTAACAGGACCCTGAGAACTTTTTTTTTTTTTTTTTTTTTTTTTTTGTGAGA                                             
# GTGTGTTTTTGTTGTTTTGACCCTTCTCACGTTGTTTTTTAGTATGCCCCATTTGTCCACAATTTTAACAGGACCCTGAGAACTTTTGTGTGTTATTTCTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTG            
#                                                       AGCCACCGTTTTAACAGGACCCTGAGAACTTTTGTGTGTTATTTCCAACTTTTAAACCTGCCATTGCCGGGCCTAGGAGGGTAGCTTTGTATACATGTCCTCCAATG
#             CCGCCCGCCTCGGCCTCCAAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCAACAGGACCCTGAGAACTTTTGTGTGTTATTTCCAACTTTTAAACCTGCCATTGCCGGGCCTAGGAGGGTAGCTTTGTATACATGTCCTCCAATG
#             CCGCCCGCCTCGGCCTCCAAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCAACCCGACCCTGAGAACTTTTGTGTGTTATTTCCAACTTTTAAACCTGCCATTGCCTGGCCTAGGAGGGTAGCTTTGTATACATTTCCTTCAATG
#
# TATCACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTGTTCCTCAAATGAGTAAAAGGCACTTCTTTACTGCTGGAATGAGGCATGTAGGTGAAAGTAAGTTGAAGTAGTCGAATAACATCTATACAGTGAGTCCTGCAAGA
# TATCACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTATTCGTAAAATGAGTAAAAGGCACTTCTGTACTGCTGGAATGAGGCATGTAGTTGAAAGTAAGTTGAAGTAGTCGAATAACATCGATCCAGTGAGTCCTGCAAGA
# TATCACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTATTCCTGAAAAGAGTAAAAGGCACTTCTGTACTGCTGGAATGAGGCATGGAGTGGAAAGGAAGTTGAAGTAGTCTGATAACAACTATCCAGGG            
# TATCACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTGTTCCTCAAATGAGTAAAAGGCACTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGGAGCGGAGGTCCTCGCT                                  
# TATCACTGTCAAAGTAAAAAACCTATTGTCCACGTCAAGGGCCAAGCTGACGTCCTGTTCCTCAAATGAGTAAAAGGCACTTCTTTACTGCTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGAATAACATCTATCCCGTGAGTCCTGCAAGA
#         ATGATCCACCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGAAATGAGTAAAAGGCACTTCTTTACTGCTGGAATGAGGCATGTAGGTGAAAGTAAGTTGAAGTAGTCGAATAACATCTATACAGTGAGTCCTGCAAGA

# Mobster input:
# #Version: 0.1.6
# #Properties file initialization, with following specified values : SAMPLENAME=test_sample MEAN_FRAGMENT_LENGTH=318 MAX_SPACING_OF_CLIPPED_READS=15 MAXIMUM_MISMATCHES_POLYA=1 MULTIPLE_SAMPLE_CALLING_STRINGENT=false DISCORDANT_CLUSTERS_MAX_OVERLAP=50 ANCHOR_SPLIT_BAM_FILE=CHM1_Illumina_all3_withoutM_sorted_markedDupl_noH_usingproperties_splitanchors.bam NEIGHBORHOOD_WINDOW_BP=200 MOBIOME_MAPPING_CMD=MosaikBuild -q (FASTQ) -st illumina -out (DAT_FILE) -quiet && MosaikAligner -in (DAT_FILE) -out (OUT_FILE) -ia /nvme/mobile_elements_2016dec/Mobster/mobiome/54_mobiles_inclHERVK.dat -hs 9 -mmp 0.1 -act 20 -j /nvme/mobile_elements_2016dec/Mobster/mobiome/54_mobiles_inclHERVK_hs9 -p 2 -annpe /nvme/mobile_elements_2016dec/Mobster/MOSAIK/2.1.26.pe.100.0065.ann -annse /nvme/mobile_elements_2016dec/Mobster/MOSAIK/2.1.26.se.100.005.ann -quiet PAIRED_END=true SD_FRAGMENT_LENGTH=103 READS_PER_CLUSTER=1 ANCHOR_BAM_FILE=CHM1_Illumina_all3_withoutM_sorted_markedDupl_noH_usingproperties_anchors.bam BAM_FROM_POTENTIALMEIFINDER=CHM1_Illumina_all3_withoutM_sorted_markedDupl_noH_usingproperties_potential.bam LENGTH_99PROCENT_OF_FRAGMENTS=901 MINIMUM_POLYA_LENGTH=9 MAPPING_TOOL=bwa IN_FILE=CHM1_Illumina_all3_withoutM_sorted_markedDupl_noH.bam MINIMUM_CLIP_LENGTH=35 MAXIMUM_CLIP_LENGTH=7 MULTIPLE_SAMPLE_CALLING=false BAM_FROM_MOBIOME_MAPPING=CHM1_Illumina_all3_withoutM_sorted_markedDupl_noH_usingproperties_mappedpotentials.bam DISCORDANT_CLUSTER_MAX_DISTANCE=600 USE_PICARD=true USE_SPLIT=true MAX_DISTANCE_OF_CLIPPED_CLUSTERS=20 MINIMUM_SUPPORTING_READS=5 MAX_OVERLAP_OF_CLIPPED_CLUSTERS=50 MOBIOME_MAPPING_TOOL=mosaik MINIMUM_AVG_QUALITY=10 REPEATMASK_FILE=/nvme/mobile_elements_2016dec/Mobster/repmask/hg19_alul1svaerv.txt MINIMUM_INITIAL_SPLIT_CLUSTER_READS=2 PICARD_COLLECT_INSERT_SIZE_METRICS_JAR=/nvme/mobile_elements_2016dec/Mobster/picard-1.73/CollectInsertSizeMetrics.jar MAX_RECORDS_IN_RAM=10000000 OUT_FILE=CHM1_Illumina_all3_withoutM_sorted_markedDupl_noH_usingproperties 
# #Creation date: Wed Dec 21 12:05:38 AEDT 2016
# Chr	Mobile Element	Insert Point	border5	border3	merged	sample	sample_counts	cluster5 length	cluster3 length	cluster5 hits	cluster3 hits	split5 hits	split3 hits	polyA5 hits	polyT5 hits	polyA3 hits	polyT3 hits	original discordant unique	original multiple	original unmapped	leftclipped max dist	rightclipped max dist	leftclipped same pos	rightclipped same pos	clipped avg qual	clipped avg length	target site duplication
# chr1	L1	10245	10225	10265	false	test_sample	test_sample=5	197	NA	5	NA	2	0	0	0	0	0	3	0	0	-1	0	0.0	1.0	33.77	40.5	unknown
# chr1	ALU	564684	564336	564704	false	test_sample	test_sample=5	NA	553	NA	5	0	0	0	0	0	0	5	0	0	-1	-1	-1	-1	-1	-1	unknown
# chr1	L1	812213	811663	812303	false	test_sample	test_sample=10	NA	281	NA	10	0	0	0	0	0	0	10	0	0	-1	-1	-1	-1	-1	-1	unknown
# chr1	ALU	3077859	3077640	3078538	false	test_sample	test_sample=8	23	NA	8	NA	0	0	0	0	0	0	8	0	0	-1	-1	-1	-1	-1	-1	unknown
# chr1	ALU	3515292	3515279	3515305	false	test_sample	test_sample=7	26	25	4	3	0	0	0	0	0	0	7	0	0	-1	-1	-1	-1	-1	-1	duplication
# chr1	ALU	7395213	7395196	7395230	false	test_sample	test_sample=13	35	32	8	5	0	0	0	0	0	0	13	0	0	-1	-1	-1	-1	-1	-1	duplication
# chr1	ALU	8580812	8580797	8580827	true	test_sample	test_sample=7	28	38	1	6	0	0	0	0	0	0	7	0	0	-1	-1	-1	-1	-1	-1	duplication
# chr1	L1	9613278	9613256	9613300	false	test_sample	test_sample=11	45	47	9	2	0	0	0	0	0	0	11	0	0	-1	-1	-1	-1	-1	-1	duplication
# chr1	ALU	12862238	12862079	12862857	false	test_sample	test_sample=7	143	NA	7	NA	0	0	0	0	0	0	7	0	0	-1	-1	-1	-1	-1	-1	unknown
# chr1	ALU	12909491	12909357	12910085	false	test_sample	test_sample=16	193	NA	16	NA	0	0	0	0	0	0	16	0	0	-1	-1	-1	-1	-1	-1	unknown
# chr1	ALU	13170390	13169765	13170555	false	test_sample	test_sample=5	NA	131	NA	5	0	0	0	0	0	0	5	0	0	-1	-1	-1	-1	-1	-1	unknown

# Output file contains same columns as input file plus following columns:
# MEI_call_quality_from_Mobster_counts
# refined_border5
# refined_border3
# refined_border5_MEI_sequence_fragment_start_pos
# refined_border3_MEI_sequence_fragment_start_pos
# refined_border5_avg_read_count_for_MEI_read_fragment
# refined_border3_avg_read_count_for_MEI_read_fragment
# refined_border5_length_of_MEI_read_fragment
# refined_border3_length_of_MEI_read_fragment
# refined_border5_percent_that_non_MEI_part_of_MEI_fragments_matches_reference
# refined_border3_percent_that_non_MEI_part_of_MEI_fragments_matches_reference
# refined_ref_seq_spanning_sequence_seen
# refined_ref_seq_spanning_sequence_seen_start_pos
# refined_ref_seq_spanning_sequence_seen_end_pos
# refined_ref_seq_spanning_sequence_seen_num_reads
# refined_border5_MEI_sequence_fragment
# refined_border3_MEI_sequence_fragment

# Mobster refers to:
# Genome Biol. 2014;15(10):488.
# Mobster: accurate detection of mobile element insertions in next generation sequencing data.
# Thung DT, de Ligt J, Vissers LE, Steehouwer M, Kroon M, de Vries P, Slagboom EP, Ye K, Veltman JA, Hehir-Kwa JY.

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import random
import subprocess
import commands
import argparse
import re
from multiprocessing import Pool


bam_reads = [] # These are the unique concatenated bam reads observed in the mapped reads bam file in the region of a mobile element insertion
bam_reads_count = [] # For each position of bam_reads, contains a 1-digit count of how many bam reads were seen having this base-pair in this position, highest value is 9
bam_reads_compare = [] # For each position of bam_reads, contains 1 or 0 to indicate whether this position matches the ref.seq. position
ref_seq_startpos = -1 # also equals bam_reads_startpos
ref_seq = '' # This is the sequence of the reference genome in the region of a mobile element insertion
read_length = 0 # This is the general length of reads in the bam file
existing_ME_regions_chrom = []
existing_ME_regions_border5 = []
existing_ME_regions_border3 = []

######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def collapse_seq( in_seq ):

	out_seq = ''
	for this_bp in in_seq:
		out_seq = out_seq + str(this_bp)
	return out_seq

######################################################
def collapse_and_strip_seq( in_seq ):

	out_seq = ''
	for this_bp in in_seq:
		out_seq = out_seq + str(this_bp)
	out_seq = out_seq.strip()
	return out_seq

######################################################
def for_qual_calc(input_count):
	output_count = input_count
	if (input_count > 40):
		output_count = 40
        return output_count

######################################################
def convert_to_integer(input_string):
	s_int = 0
	if ((input_string != 'NA') and (input_string != '')):
		s_int = int(input_string)
        return s_int

######################################################
def parse_out_any_warnings( command_output, command ):

	outlines = []

	all_outlines = command_output.split("\n")
	for outline in all_outlines:
		is_a_warning_line = False
		if (len(outline) >= 7):
			if (outline[0:7] == 'Warning'):
				is_a_warning_line = True
		if (is_a_warning_line == False):
			outlines.append( outline )

	return outlines

######################################################
def calculate_average_counts( count_array ):

	sum_of_counts = 0
	number_of_bp = 0
	for this_bp in count_array:
		if ((this_bp != ' ') and (this_bp != '0')):
			sum_of_counts = sum_of_counts + int(this_bp)
			number_of_bp = number_of_bp + 1
	average_count = 0
	if (number_of_bp > 0):
		average_count = float(sum_of_counts) / float(number_of_bp)
	return average_count

######################################################
def calculate_mismatch_region_length_from_the_left_end( aligned_match_sequence, minimum_mismatch_percent ):

	region_length = 0
	region_length_1by1 = 0
	region_length_byWindow = 0
	match_sequence = aligned_match_sequence.strip()

	keep_looking = True
	bp_upto = 1
	if (bp_upto > len(match_sequence)):
		keep_looking = False
	while (keep_looking):
		count_mismatches = match_sequence.count( '0', 0, bp_upto )
		percent_mismatches = count_mismatches * 100 / bp_upto
		if (percent_mismatches >= minimum_mismatch_percent):
			region_length_1by1 = bp_upto
			bp_upto = bp_upto + 1
			if (bp_upto > len(match_sequence)):
				keep_looking = False
		else:
			keep_looking = False

	window_length = 10
	keep_looking = True
	bp_upto = window_length
	if (bp_upto > len(match_sequence)):
		keep_looking = False
	while (keep_looking):
		count_mismatches = match_sequence.count( '0', 0, bp_upto )
		percent_mismatches = count_mismatches * 100 / bp_upto
		if (percent_mismatches >= minimum_mismatch_percent):
			region_length_byWindow = bp_upto
			bp_upto = bp_upto + window_length
			if (bp_upto > len(match_sequence)):
				keep_looking = False
		else:
			keep_looking = False

	region_length = max( region_length_1by1, region_length_byWindow )

	return region_length

######################################################
def calculate_mismatch_region_length_from_the_right_end( aligned_match_sequence, minimum_mismatch_percent ):

	region_length = 0
	region_length_1by1 = 0
	region_length_byWindow = 0
	match_sequence = aligned_match_sequence.strip()

	keep_looking = True
	bp_upto = len(match_sequence) - 1
	if (bp_upto <= 0):
		keep_looking = False
	while (keep_looking):
		count_mismatches = match_sequence.count( '0', bp_upto, len(match_sequence) )
		percent_mismatches = 0
		denom = len(match_sequence) - bp_upto
		if (denom > 0):
			percent_mismatches = count_mismatches * 100 / denom
		if (percent_mismatches >= minimum_mismatch_percent):
			region_length_1by1 = len(match_sequence) - bp_upto
			bp_upto = bp_upto - 1
			if (bp_upto <= 0):
				keep_looking = False
		else:
			keep_looking = False

	window_length = 10
	keep_looking = True
	bp_upto = len(match_sequence) - window_length
	if (bp_upto <= 0):
		keep_looking = False
	while (keep_looking):
		count_mismatches = match_sequence.count( '0', bp_upto, len(match_sequence) )
		percent_mismatches = 0
		denom = len(match_sequence) - bp_upto
		if (denom > 0):
			percent_mismatches = count_mismatches * 100 / denom
		if (percent_mismatches >= minimum_mismatch_percent):
			region_length_byWindow = len(match_sequence) - bp_upto
			bp_upto = bp_upto - window_length
			if (bp_upto <= 0):
				keep_looking = False
		else:
			keep_looking = False

	region_length = max( region_length_1by1, region_length_byWindow )

	return region_length

######################################################
def find_3prime_end_of_mobile_element_insertion(): # looking for inhouse_border5

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq
	return_result = ''
	found_the_mobile_element_end = False
	continue_looking_at_reads = True
	read_upto = 0
	if (read_upto >= len(bam_reads)):
		continue_looking_at_reads = False
	sliding_window = 10
	minimum_match_percent = 90
	minimum_mismatch_percent = 50
	minimum_match_amount = int(sliding_window * minimum_match_percent / 100)
	minimum_mismatch_amount = int(sliding_window * minimum_mismatch_percent / 100)
	best_pos = -1
	best_percent_match_on_the_right = -1
	best_seq = ''
	best_seq_stripped = ''
	best_seq_start_pos = -1
	best_average_counts_on_the_left = -1
	best_mismatch_region_length_on_left = -1

	while (continue_looking_at_reads):

		this_read = bam_reads[read_upto]
		this_read_count = bam_reads_count[read_upto]
		this_read_compare = bam_reads_compare[read_upto]
		this_read_compare_collapsed = collapse_seq(this_read_compare)

		continue_looking_at_positions = True
		pos_upto = len(this_read) - 1
		found_end_of_this_seq = False
		while (found_end_of_this_seq == False):
			if ((this_read[pos_upto] != ' ') and (this_read[pos_upto] != 'N')):
				found_end_of_this_seq = True
			else:
				pos_upto = pos_upto - 1
				if (pos_upto < 0):
					found_end_of_this_seq = True
		pos_upto = pos_upto - sliding_window
		if (pos_upto < 0):
			continue_looking_at_positions = False
		found_matching_region = False
		found_mismatching_region = False
		rightmost_matching_region = ''
		leftmost_mismatching_region = ''
		# Start from the right, look for the left-most end of the matching region, having mismatching region to the left of it. That is the mobile element insert point.
		while (continue_looking_at_positions):

			this_read_pos_bp = this_read[pos_upto]
			this_ref_pos_bp = ref_seq[pos_upto]

			# if ( (found_matching_region == False) or ((found_matching_region == True) and (found_mismatching_region == False)) ):
			if (found_mismatching_region == False):
				if (this_read_pos_bp == this_ref_pos_bp):
					pos1 = pos_upto
					pos2 = pos_upto + sliding_window
					num_match = this_read_compare_collapsed.count( '1', pos1, pos2 ) # compare the sliding window to the right of and including pos_upto
					if (num_match >= minimum_match_amount):
						found_matching_region = True
						rightmost_matching_region = pos_upto

			if (found_matching_region):
				next_pos_upto = pos_upto - 1
				if (next_pos_upto >= 0):
					next_read_pos_bp = this_read[next_pos_upto]
					next_ref_pos_bp = ref_seq[next_pos_upto]
					if ((next_read_pos_bp != next_ref_pos_bp) and (next_read_pos_bp != ' ') and (next_ref_pos_bp != ' ')):
						pos1 = pos_upto - sliding_window
						pos2 = pos_upto
						if (pos1 >= 0):
							num_mismatch = this_read_compare_collapsed.count( '0', pos1, pos2 ) # compare the sliding window to the left of pos_upto ending at the position before pos_upto
							if (num_mismatch >= minimum_mismatch_amount):
								found_mismatching_region = True
								leftmost_mismatching_region = next_pos_upto

			if ((found_mismatching_region) and (found_matching_region)):

				this_seq = collapse_seq(this_read)
				this_seq_stripped = this_seq.strip()
				this_seq_start_pos = this_seq.find( this_seq_stripped )
				pos1 = pos_upto
				pos2 = pos_upto + 50 # don't use pos2 = len(this_read) because that makes this_percent_match_on_the_right dependent on concatenated read length, not on how well it matches to the right of this position 

				# how much is this insertion on the right of this position supported by a lot of reads showing this alternate sequence on the right
				this_average_counts_on_the_left = calculate_average_counts( this_read_count[ 0 : pos_upto ] )
				percent_increase_in_average_counts_on_the_left = -1
				this_mismatch_region_length_on_left = -1
				percent_increase_in_mismatch_region_on_left = -1
				this_percent_match_on_the_right = -1

				# Various criteria to consider when deciding if this sequence represents a better representation of this mobile-element-insertion (MEI) than the previous sequence we found to represent this MEI.
				# If this sequence has more than 1 supporting bam read (this_average_counts_on_the_left/right) ==> weak.
				# If this sequence has more than 1 supporting bam read (this_average_counts_on_the_left/right) ==> good.
				# If this sequence has quite a few supporting bam reads (this_average_counts_on_the_left/right) ==> very good.
				# If this sequence has less than 10 mobile-element-insertion base-pairs (this_mismatch_region_on_the_left/right) ==> weak.
				# If this sequence has more than 10 mobile-element-insertion base-pairs (this_mismatch_region_on_the_left/right) ==> good.

				choose_this_sequence = False # there may not be a best_pos yet
				if ((this_average_counts_on_the_left > 1) and (best_average_counts_on_the_left <= 1)):
					choose_this_sequence = True
				elif ((this_average_counts_on_the_left <= 1) and (best_average_counts_on_the_left > 1)):
					choose_this_sequence = False
				elif ((this_average_counts_on_the_left > 1) and (best_average_counts_on_the_left > 1)):

					# if need to choose between this read and a previously found read, then
					# if both have > 1 supporting read, and both have mismatch length >= 10, then choose based on how well the matching part of the read matches the ref.seq. 
					# if both have > 1 supporting read, but not both have mismatch length >= 10, then choose the one that has mismatch length >= 10
					# if both have > 1 supporting read, and both have mismatch length < 10, then choose based on how well the matching part of the read matches the ref.seq. 

					# how long is this mismatching region on the left of this position
					this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )
					if ((this_mismatch_region_length_on_left >= 10) and (best_mismatch_region_length_on_left >= 10)):
						# how well does this read match the ref.seq. on the right before this pos_upto after which it doesn't match
						this_percent_match_on_the_right = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
						if (this_percent_match_on_the_right > best_percent_match_on_the_right):
							choose_this_sequence = True
					elif ((this_mismatch_region_length_on_left >= 10) and (best_mismatch_region_length_on_left < 10)):
						choose_this_sequence = True
					elif ((this_mismatch_region_length_on_left < 10) and (best_mismatch_region_length_on_left < 10)):
						# how well does this read match the ref.seq. on the right before this pos_upto after which it doesn't match
						this_percent_match_on_the_right = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
						if (this_percent_match_on_the_right > best_percent_match_on_the_right):
							choose_this_sequence = True

				else: # ((this_average_counts_on_the_left <= 1) and (best_average_counts_on_the_left <= 1)), there may be no best_pos yet
					# both sequences are poor representations of this MEI when considering only number of supporting reads, both only have 1 supporting read

					# if need to choose between this read and a previously found read, then
					# if both have <= 1 supporting read, and both have mismatch length >= 10, then choose based on how well the matching part of the read matches the ref.seq. 
					# if both have <= 1 supporting read, but not both have mismatch length >= 10, then choose the one that has mismatch length >= 10
					# if both have <= 1 supporting read, and both have mismatch length < 10, then choose based on how well the matching part of the read matches the ref.seq. 

					if (best_pos == -1):
						choose_this_sequence = True
					else:
						# how long is this mismatching region on the left of this position
						this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )
						if ((this_mismatch_region_length_on_left >= 10) and (best_mismatch_region_length_on_left >= 10)):
							# how well does this read match the ref.seq. on the right before this pos_upto after which it doesn't match
							this_percent_match_on_the_right = 0
							denom = pos2 - pos1 + 1
							if (denom > 0):
								this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
							if (this_percent_match_on_the_right > best_percent_match_on_the_right):
								choose_this_sequence = True
						elif ((this_mismatch_region_length_on_left >= 10) and (best_mismatch_region_length_on_left < 10)):
							choose_this_sequence = True
						elif ((this_mismatch_region_length_on_left < 10) and (best_mismatch_region_length_on_left < 10)):
							# how well does this read match the ref.seq. on the right before this pos_upto after which it doesn't match
							this_percent_match_on_the_right = 0
							denom = pos2 - pos1 + 1
							if (denom > 0):
								this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
							if (this_percent_match_on_the_right > best_percent_match_on_the_right):
								choose_this_sequence = True

				if (choose_this_sequence):
					if (this_mismatch_region_length_on_left == -1):
						this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )
					if (this_percent_match_on_the_right == -1):
						this_percent_match_on_the_right = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
					best_pos = pos_upto
					best_seq = this_seq
					best_seq_stripped = this_seq_stripped
					best_seq_start_pos = this_seq_start_pos
					best_percent_match_on_the_right = this_percent_match_on_the_right
					best_average_counts_on_the_left = this_average_counts_on_the_left
					best_mismatch_region_length_on_left = this_mismatch_region_length_on_left

				# We found a possible MEI position in this sequence, whether or not we decided to use it.
				# Continue looking along this sequence in case there are even better possible MEI positions in this same sequence.
				# (If we don't do this, then when this sequence contains a strong MEI further along, 
				# the current weak MEI will be replaced by some other sequence's weak MEI, and we'll miss the strong MEI.)
				pos_upto = pos_upto - 1
				found_matching_region = False
				found_mismatching_region = False
				rightmost_matching_region = ''
				leftmost_mismatching_region = ''
				if (pos_upto < 0):
					continue_looking_at_positions = False

			else:
				pos_upto = pos_upto - 1
				if (pos_upto < 0):
					continue_looking_at_positions = False

		read_upto = read_upto + 1
		if (read_upto >= len(bam_reads)):
			continue_looking_at_reads = False

	if (best_pos != -1):
		return_result = str(best_pos) + ',' + str(best_seq_start_pos) + ',' + str(best_average_counts_on_the_left) + ',' + str(best_mismatch_region_length_on_left) + ',' + str(best_percent_match_on_the_right) + ',' + best_seq_stripped
        return return_result

######################################################
def find_3prime_end_of_mobile_element_insertion_got_5prime_already( the_5prime_end ): # looking for inhouse_border5

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq
	return_result = ''
	found_the_mobile_element_end = False
	continue_looking_at_reads = True
	read_upto = 0
	if (read_upto >= len(bam_reads)):
		continue_looking_at_reads = False
	sliding_window = 10
	minimum_match_percent = 90
	minimum_mismatch_percent = 50
	minimum_match_amount = int(sliding_window * minimum_match_percent / 100)
	minimum_mismatch_amount = int(sliding_window * minimum_mismatch_percent / 100)
	best_pos = -1
	best_percent_match_on_the_right = -1
	best_seq = ''
	best_seq_stripped = ''
	best_seq_start_pos = -1
	best_average_counts_on_the_left = -1
	best_mismatch_region_length_on_left = -1

	while (continue_looking_at_reads):

		this_read = bam_reads[read_upto]
		this_read_count = bam_reads_count[read_upto]
		this_read_compare = bam_reads_compare[read_upto]
		this_read_compare_collapsed = collapse_seq(this_read_compare)

		continue_looking_at_positions = True
		pos_upto = len(this_read) - 1
		found_end_of_this_seq = False
		while (found_end_of_this_seq == False):
			if ((this_read[pos_upto] != ' ') and (this_read[pos_upto] != 'N')):
				found_end_of_this_seq = True
			else:
				pos_upto = pos_upto - 1
				if (pos_upto < 0):
					found_end_of_this_seq = True
		pos_upto = pos_upto - sliding_window
		if (pos_upto < 0):
			continue_looking_at_positions = False
		found_matching_region = False
		found_mismatching_region = False
		rightmost_matching_region = ''
		leftmost_mismatching_region = ''
		# Start from the right, look for the left-most end of the matching region, having mismatching region to the left of it. That is the mobile element insert point.
		while (continue_looking_at_positions):

			this_read_pos_bp = this_read[pos_upto]
			this_ref_pos_bp = ref_seq[pos_upto]

			# if ( (found_matching_region == False) or ((found_matching_region == True) and (found_mismatching_region == False)) ):
			if (found_mismatching_region == False):
				if (this_read_pos_bp == this_ref_pos_bp):
					pos1 = pos_upto
					pos2 = pos_upto + sliding_window
					num_match = this_read_compare_collapsed.count( '1', pos1, pos2 ) # compare the sliding window to the right of and including pos_upto
					if (num_match >= minimum_match_amount):
						found_matching_region = True
						rightmost_matching_region = pos_upto

			if (found_matching_region):
				next_pos_upto = pos_upto - 1
				if (next_pos_upto >= 0):
					next_read_pos_bp = this_read[next_pos_upto]
					next_ref_pos_bp = ref_seq[next_pos_upto]
					if ((next_read_pos_bp != next_ref_pos_bp) and (next_read_pos_bp != ' ') and (next_ref_pos_bp != ' ')):
						pos1 = pos_upto - sliding_window
						pos2 = pos_upto
						if (pos1 >= 0):
							num_mismatch = this_read_compare_collapsed.count( '0', pos1, pos2 ) # compare the sliding window to the left of pos_upto ending at the position before pos_upto			
							if (num_mismatch >= minimum_mismatch_amount):
								found_mismatching_region = True
								leftmost_mismatching_region = next_pos_upto

			if ((found_mismatching_region) and (found_matching_region)):

				this_seq = collapse_seq(this_read)
				this_seq_stripped = this_seq.strip()
				this_seq_start_pos = this_seq.find( this_seq_stripped )
				pos1 = pos_upto
				pos2 = pos_upto + 50 # don't use pos2 = len(this_read) because that makes this_percent_match_on_the_right dependent on concatenated read length, not on how well it matches to the right of this position 

				this_average_counts_on_the_left = -1
				this_mismatch_region_length_on_left = -1
				this_percent_match_on_the_right = -1

				# If there is more than one 3-prime end in the vicinity, 
				# then go for the one whose border5 position is before the border3 position of the 5-prime end found previously.
				# pos_upto and the_5prime_end are zero-based index

				choose_this_sequence = False
				if (best_pos < the_5prime_end):
					if (pos_upto > the_5prime_end): # (pos_upto < the_5prime_end) and (best_pos > the_5prime_end) so this pos_upto is clearly better than best_pos
						choose_this_sequence = True
					else: # pos_upto is also > the_5prime_end, as is best_pos. We would prefer a best_pos that is before the_5prime_end
						this_diff = abs(the_5prime_end - pos_upto)
						best_diff = abs(the_5prime_end - best_pos)
						if (this_diff < best_diff): # pick the border5 3prime-end that is closest to border3 5prime-end, unless the distances are almost the same (within 10% of each other)
							percent_increase_in_distance_from_previously_found_border = 0
							if (best_diff > 0):
								percent_increase_in_distance_from_previously_found_border = abs(this_diff - best_diff) * 100 / best_diff
							if (percent_increase_in_distance_from_previously_found_border > 10):
								choose_this_sequence = True
							else: # (percent_increase_in_distance_from_previously_found_border < 10)
								this_average_counts_on_the_left = calculate_average_counts( this_read_count[ 0 : pos_upto ] )
								# if the difference in distance is only 10% or less, then choose based on number of supporting reads or else length of mismatch region
								if ((this_average_counts_on_the_left > 1) and (best_average_counts_on_the_left <= 1)):
									choose_this_sequence = True
								elif ((this_average_counts_on_the_left <= 1) and (best_average_counts_on_the_left > 1)):
									choose_this_sequence = False
								else: # they both have > 1 supporting read, or they both have only 1 supporting read, so choose the longest mismatch region
									this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )
									if (this_mismatch_region_length_on_left > best_mismatch_region_length_on_left):
										choose_this_sequence = True
				else: # best_pos < the_5prime_end. We want a best_pos < the_5prime_end.
					if (pos_upto < the_5prime_end):
						this_diff = abs(pos_upto - the_5prime_end)
						best_diff = abs(best_pos - the_5prime_end)
						if ((this_diff <= 100) and (best_diff <= 100)): # both are close to the previously found border, so both are very good, choose based on mismatch_region, then supporting read counts, then matching ref.seq.

							# how long is this mismatching region on the left of this position
							this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )

							if ((this_mismatch_region_length_on_left >= 10) and (best_mismatch_region_length_on_left < 10)):
								choose_this_sequence = True
							elif ((this_mismatch_region_length_on_left < 10) and (best_mismatch_region_length_on_left >= 10)):
								choose_this_sequence = False
							elif ((this_mismatch_region_length_on_left >= 10) and (best_mismatch_region_length_on_left >= 10)): 
								# both sequences are good representations of this MEI when considering only length of mismatch_region_length_on_left

								this_average_counts_on_the_left = calculate_average_counts( this_read_count[ 0 : pos_upto ] )
								# if the difference in distance is only 10% or less, then choose based on number of supporting reads or else length of mismatch region
								if ((this_average_counts_on_the_left > 1) and (best_average_counts_on_the_left <= 1)):
									choose_this_sequence = True
								elif ((this_average_counts_on_the_left <= 1) and (best_average_counts_on_the_left > 1)):
									choose_this_sequence = False
								else: # they both have > 1 supporting read, or they both have only 1 supporting read, so choose using another criteria

									# how much increase in this mismatching region on the left of this position
									percent_increase_in_mismatch_region_on_left = 0
									if (best_pos == -1):
										percent_increase_in_mismatch_region_on_left = 100
									else:
										if ((best_mismatch_region_length_on_left == 0) and (this_mismatch_region_length_on_left > 0)):
											percent_increase_in_mismatch_region_on_left = 100
										else:
											if (this_mismatch_region_length_on_left > best_mismatch_region_length_on_left):
												percent_increase_in_mismatch_region_on_left = ((this_mismatch_region_length_on_left - best_mismatch_region_length_on_left) * 100) / best_mismatch_region_length_on_left

									if (percent_increase_in_mismatch_region_on_left >= 10):
										choose_this_sequence = True
									else: # (percent_increase_in_mismatch_region_on_left < 10%) so we consider the 2 values to be roughly equal

										if (this_average_counts_on_the_left > best_average_counts_on_the_left):
											choose_this_sequence = True
										elif (this_mismatch_region_length_on_left > best_mismatch_region_length_on_left):
											choose_this_sequence = True
										else:

											# how well does this read match the ref.seq. on the right before this pos_upto after which it doesn't match
											this_percent_match_on_the_right = 0
											denom = pos2 - pos1 + 1
											if (denom > 0):
												this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom

											if (this_percent_match_on_the_right > best_percent_match_on_the_right):
												choose_this_sequence = True
											elif (best_average_counts_on_the_left > this_average_counts_on_the_left):
												choose_this_sequence = False
											elif (best_mismatch_region_length_on_left > this_mismatch_region_length_on_left):
												choose_this_sequence = False
											elif (best_percent_match_on_the_right > this_percent_match_on_the_right):
												choose_this_sequence = False
											else: # the 2 sequences are just as good or bad as each other as representing this MEI
												choose_this_sequence = False

							else: # ((this_mismatch_region_length_on_left < 10 base-pairs) and (best_mismatch_region_length_on_left < 10 base-pairs))
								# both sequences are poor representations of this MEI when considering only length of mismatch_region_length_on_left
								# however, they both have more than one supporting read

								if (this_average_counts_on_the_left > best_average_counts_on_the_left):
									choose_this_sequence = True
								elif (this_mismatch_region_length_on_left > best_mismatch_region_length_on_left):
									choose_this_sequence = True
								else:

									# how well does this read match the ref.seq. on the right before this pos_upto after which it doesn't match
									this_percent_match_on_the_right = 0
									denom = pos2 - pos1 + 1
									if (denom > 0):
										this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom

									if (this_percent_match_on_the_right > best_percent_match_on_the_right):
										choose_this_sequence = True
									elif (best_average_counts_on_the_left > this_average_counts_on_the_left):
										choose_this_sequence = False
									elif (best_mismatch_region_length_on_left > this_mismatch_region_length_on_left):
										choose_this_sequence = False
									elif (best_percent_match_on_the_right > this_percent_match_on_the_right):
										choose_this_sequence = False
									else: # the 2 sequences are just as good or bad as each other as representing this MEI
										choose_this_sequence = False

						else:
							if (this_diff < best_diff): # pick the border5 3prime-end that is closest to border3 5prime-end, unless the distances are almost the same (within 10% of each other)
								percent_increase_in_distance_from_previously_found_border = 0
								if (best_diff > 0):
									percent_increase_in_distance_from_previously_found_border = abs(this_diff - best_diff) * 100 / best_diff
								if (percent_increase_in_distance_from_previously_found_border > 10):
									choose_this_sequence = True
								else: # (percent_increase_in_distance_from_previously_found_border < 10)

									this_average_counts_on_the_left = calculate_average_counts( this_read_count[ 0 : pos_upto ] )
									# how much increase in the number of reads showing this alternate sequence on the right
									percent_increase_in_average_counts_on_the_left = 0
									if (best_pos == -1):
										percent_increase_in_average_counts_on_the_left = 100
									else:
										if ((best_average_counts_on_the_left == 0) and (this_average_counts_on_the_left > 0)):
											percent_increase_in_average_counts_on_the_left = 100
										else:
											if (this_average_counts_on_the_left > best_average_counts_on_the_left):
												percent_increase_in_average_counts_on_the_left = 0
												if (best_average_counts_on_the_left > 0):
													percent_increase_in_average_counts_on_the_left = ((this_average_counts_on_the_left - best_average_counts_on_the_left) * 100) / best_average_counts_on_the_left

									if (percent_increase_in_average_counts_on_the_left >= 10):
										choose_this_sequence = True
									else: # (percent_increase_in_average_counts_on_the_left < 10%) so we consider the 2 values to be roughly equal

										# if the difference in distance is only 10% or less, then choose based on number of supporting reads or else length of mismatch region
										if ((this_average_counts_on_the_left > 1) and (best_average_counts_on_the_left <= 1)):
											choose_this_sequence = True
										elif ((this_average_counts_on_the_left <= 1) and (best_average_counts_on_the_left > 1)):
											choose_this_sequence = False
										else: # they both have > 1 supporting read, or they both have only 1 supporting read, so choose the longest mismatch region
											this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )
											if (this_mismatch_region_length_on_left > best_mismatch_region_length_on_left):
												choose_this_sequence = True
					# else (pos_upto >= the_5prime_end) and (best_pos < the_5prime_end), so keep the current best_pos, it is better than this pos_upto

				if (choose_this_sequence):
					if (this_average_counts_on_the_left == -1):
						this_average_counts_on_the_left = calculate_average_counts( this_read_count[ 0 : pos_upto ] )
					if (this_mismatch_region_length_on_left == -1):
						this_mismatch_region_length_on_left = calculate_mismatch_region_length_from_the_right_end( this_read_compare_collapsed[ 0 : pos_upto ], minimum_mismatch_percent )
					if (this_percent_match_on_the_right == -1):
						this_percent_match_on_the_right = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_right = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
					best_pos = pos_upto
					best_seq = this_seq
					best_seq_stripped = this_seq_stripped
					best_seq_start_pos = this_seq_start_pos
					best_percent_match_on_the_right = this_percent_match_on_the_right
					best_average_counts_on_the_left = this_average_counts_on_the_left
					best_mismatch_region_length_on_left = this_mismatch_region_length_on_left

				# We found a possible MEI position in this sequence, whether or not we decided to use it.
				# Continue looking along this sequence in case there are even better possible MEI positions in this same sequence.
				# (If we don't do this, then when this sequence contains a strong MEI further along, 
				# the current weak MEI will be replaced by some other sequence's weak MEI, and we'll miss the strong MEI.)
				pos_upto = pos_upto - 1
				found_matching_region = False
				found_mismatching_region = False
				rightmost_matching_region = ''
				leftmost_mismatching_region = ''
				if (pos_upto < 0):
					continue_looking_at_positions = False

			else:
				pos_upto = pos_upto - 1
				if (pos_upto < 0):
					continue_looking_at_positions = False

		read_upto = read_upto + 1
		if (read_upto >= len(bam_reads)):
			continue_looking_at_reads = False

	if (best_pos != -1):
		return_result = str(best_pos) + ',' + str(best_seq_start_pos) + ',' + str(best_average_counts_on_the_left) + ',' + str(best_mismatch_region_length_on_left) + ',' + str(best_percent_match_on_the_right) + ',' + best_seq_stripped
        return return_result

######################################################
def find_5prime_end_of_mobile_element_insertion(): # looking for inhouse_border3

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq
	return_result = ''
	found_the_mobile_element_end = False
	continue_looking_at_reads = True
	read_upto = 0
	if (read_upto >= len(bam_reads)):
		continue_looking_at_reads = False
	sliding_window = 10
	minimum_match_percent = 90
	minimum_mismatch_percent = 50
	minimum_match_amount = int(sliding_window * minimum_match_percent / 100)
	minimum_mismatch_amount = int(sliding_window * minimum_mismatch_percent / 100)
	best_pos = -1
	best_percent_match_on_the_left = -1
	best_seq = ''
	best_seq_stripped = ''
	best_seq_start_pos = -1
	best_average_counts_on_the_right = -1
	best_mismatch_region_length_on_right = -1

	while (continue_looking_at_reads):

		this_read = bam_reads[read_upto]
		this_read_count = bam_reads_count[read_upto]
		this_read_compare = bam_reads_compare[read_upto]
		this_read_compare_collapsed = collapse_seq(this_read_compare)

		continue_looking_at_positions = True
		pos_upto = sliding_window - 1
		if (pos_upto >= len(this_read)):
			continue_looking_at_positions = False
		found_matching_region = False
		found_mismatching_region = False
		rightmost_matching_region = ''
		leftmost_mismatching_region =''
		# Start from the left, look for the right-most end of the matching region, having mismatching region to the right of it. That is the mobile element insert point.
		while (continue_looking_at_positions):

			this_read_pos_bp = this_read[pos_upto]
			this_ref_pos_bp = ref_seq[pos_upto]

			# if ( (found_matching_region == False) or ((found_matching_region == True) and (found_mismatching_region == False)) ):
			if (found_mismatching_region == False):
				if (this_read_pos_bp == this_ref_pos_bp):
					pos1 = pos_upto - sliding_window + 1
					pos2 = pos_upto + 1
					num_match = this_read_compare_collapsed.count( '1', pos1, pos2 ) # compare the sliding window to the left of and including up to pos_upto
					if (num_match >= minimum_match_amount):
						found_matching_region = True
						rightmost_matching_region = pos_upto

			if (found_matching_region):
				next_pos_upto = pos_upto + 1
				if (next_pos_upto < len(this_read)):
					next_read_pos_bp = this_read[next_pos_upto]
					next_ref_pos_bp = ref_seq[next_pos_upto]
					if ((next_read_pos_bp != next_ref_pos_bp) and (next_read_pos_bp != ' ') and (next_ref_pos_bp != ' ')):
						pos1 = pos_upto + 1
						pos2 = pos_upto + sliding_window
						if (pos2 < len(this_read)):
							num_mismatch = this_read_compare_collapsed.count( '0', pos1, pos2 ) # compare the sliding window to the right of pos_upto starting at the position after pos_upto			
							if (num_mismatch >= minimum_mismatch_amount):
								found_mismatching_region = True
								leftmost_mismatching_region = next_pos_upto

			if ((found_mismatching_region) and (found_matching_region)):

				this_seq = collapse_seq(this_read)
				this_seq_stripped = this_seq.strip()
				this_seq_start_pos = this_seq.find( this_seq_stripped )
				pos2 = pos_upto + 1
				pos1 = pos_upto - 50 # don't use pos1 = 0 because that makes this_percent_match_on_the_right dependent on concatenated read length, not on how well it matches to the left of this position 

				# how much is this insertion on the left of this position supported by a lot of reads showing this alternate sequence on the left
				this_average_counts_on_the_right = calculate_average_counts( this_read_count[ 0 : pos_upto ] )
				percent_increase_in_average_counts_on_the_right = -1
				this_mismatch_region_length_on_right = -1
				percent_increase_in_mismatch_region_on_right = -1
				this_percent_match_on_the_left = -1

				# Various criteria to consider when deciding if this sequence represents a better representation of this mobile-element-insertion (MEI) than the previous sequence we found to represent this MEI.
				# If this sequence has more than 1 supporting bam read (this_average_counts_on_the_left/right) ==> weak.
				# If this sequence has more than 1 supporting bam read (this_average_counts_on_the_left/right) ==> good.
				# If this sequence has quite a few supporting bam reads (this_average_counts_on_the_left/right) ==> very good.
				# If this sequence has less than 10 mobile-element-insertion base-pairs (this_mismatch_region_on_the_left/right) ==> weak.
				# If this sequence has more than 10 mobile-element-insertion base-pairs (this_mismatch_region_on_the_left/right) ==> good.

				choose_this_sequence = False # there may not be a best_pos yet
				if ((this_average_counts_on_the_right > 1) and (best_average_counts_on_the_right <= 1)):
					choose_this_sequence = True
				elif ((this_average_counts_on_the_right <= 1) and (best_average_counts_on_the_right > 1)):
					choose_this_sequence = False
				elif ((this_average_counts_on_the_right > 1) and (best_average_counts_on_the_right > 1)):

					# if need to choose between this read and a previously found read, then
					# if both have > 1 supporting read, and both have mismatch length >= 10, then choose based on how well the matching part of the read matches the ref.seq. 
					# if both have > 1 supporting read, but not both have mismatch length >= 10, then choose the one that has mismatch length >= 10
					# if both have > 1 supporting read, and both have mismatch length < 10, then choose based on how well the matching part of the read matches the ref.seq. 

					# how long is this mismatching region on the right of this position
					this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )
					if ((this_mismatch_region_length_on_right >= 10) and (best_mismatch_region_length_on_right >= 10)):
						# how well does this read match the ref.seq. on the left before this pos_upto after which it doesn't match
						this_percent_match_on_the_left = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
						if (this_percent_match_on_the_left > best_percent_match_on_the_left):
							choose_this_sequence = True
					elif ((this_mismatch_region_length_on_right >= 10) and (best_mismatch_region_length_on_right < 10)):
						choose_this_sequence = True
					elif ((this_mismatch_region_length_on_right < 10) and (best_mismatch_region_length_on_right < 10)):
						# how well does this read match the ref.seq. on the left before this pos_upto after which it doesn't match
						this_percent_match_on_the_left = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
						if (this_percent_match_on_the_left > best_percent_match_on_the_left):
							choose_this_sequence = True

				else: # ((this_average_counts_on_the_right <= 1) and (best_average_counts_on_the_right <= 1)), there may be no best_pos yet
					# both sequences are poor representations of this MEI when considering only number of supporting reads, both only have 1 supporting read

					if (best_pos == -1):
						choose_this_sequence = True
					else:

						# if need to choose between this read and a previously found read, then
						# if both have <= 1 supporting read, and both have mismatch length >= 10, then choose based on how well the matching part of the read matches the ref.seq. 
						# if both have <= 1 supporting read, but not both have mismatch length >= 10, then choose the one that has mismatch length >= 10
						# if both have <= 1 supporting read, and both have mismatch length < 10, then choose based on how well the matching part of the read matches the ref.seq. 

						# how long is this mismatching region on the right of this position
						this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )
						if ((this_mismatch_region_length_on_right >= 10) and (best_mismatch_region_length_on_right >= 10)):
							# how well does this read match the ref.seq. on the left before this pos_upto after which it doesn't match
							this_percent_match_on_the_left = 0
							denom = pos2 - pos1 + 1
							if (denom > 0):
								this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
							if (this_percent_match_on_the_left > best_percent_match_on_the_left):
								choose_this_sequence = True
						elif ((this_mismatch_region_length_on_right >= 10) and (best_mismatch_region_length_on_right < 10)):
							choose_this_sequence = True
						elif ((this_mismatch_region_length_on_right < 10) and (best_mismatch_region_length_on_right < 10)):
							# how well does this read match the ref.seq. on the left before this pos_upto after which it doesn't match
							this_percent_match_on_the_left = 0
							denom = pos2 - pos1 + 1
							if (denom > 0):
								this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
							if (this_percent_match_on_the_left > best_percent_match_on_the_left):
								choose_this_sequence = True

				if (choose_this_sequence):
					if (this_mismatch_region_length_on_right == -1):
						this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )
					if (this_percent_match_on_the_left == -1):
						this_percent_match_on_the_left = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
					best_pos = pos_upto
					best_seq = this_seq
					best_seq_stripped = this_seq_stripped
					best_seq_start_pos = this_seq_start_pos
					best_percent_match_on_the_left = this_percent_match_on_the_left
					best_average_counts_on_the_right = this_average_counts_on_the_right
					best_mismatch_region_length_on_right = this_mismatch_region_length_on_right

				# We found a possible MEI position in this sequence, whether or not we decided to use it.
				# Continue looking along this sequence in case there are even better possible MEI positions in this same sequence.
				# (If we don't do this, then when this sequence contains a strong MEI further along, 
				# the current weak MEI will be replaced by some other sequence's weak MEI, and we'll miss the strong MEI.)
				pos_upto = pos_upto + 1
				found_matching_region = False
				found_mismatching_region = False
				rightmost_matching_region = ''
				leftmost_mismatching_region =''
				if (pos_upto >= len(this_read)):
					continue_looking_at_positions = False

			else:
				pos_upto = pos_upto + 1
				if (pos_upto >= len(this_read)):
					continue_looking_at_positions = False

		read_upto = read_upto + 1
		if (read_upto >= len(bam_reads)):
			continue_looking_at_reads = False

	if (best_pos != -1):
		return_result = str(best_pos) + ',' + str(best_seq_start_pos) + ',' + str(best_average_counts_on_the_right) + ',' + str(best_mismatch_region_length_on_right) + ',' + str(best_percent_match_on_the_left) + ',' + best_seq_stripped
        return return_result

######################################################
def find_5prime_end_of_mobile_element_insertion_got_3prime_already( the_3prime_end ): # looking for inhouse_border3

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq
	return_result = ''
	found_the_mobile_element_end = False
	continue_looking_at_reads = True
	read_upto = 0
	if (read_upto >= len(bam_reads)):
		continue_looking_at_reads = False
	sliding_window = 10
	minimum_match_percent = 90
	minimum_mismatch_percent = 50
	minimum_match_amount = int(sliding_window * minimum_match_percent / 100)
	minimum_mismatch_amount = int(sliding_window * minimum_mismatch_percent / 100)
	best_pos = -1
	best_percent_match_on_the_left = -1
	best_seq = ''
	best_seq_stripped = ''
	best_seq_start_pos = -1
	best_average_counts_on_the_right = -1
	best_mismatch_region_length_on_right = -1

	while (continue_looking_at_reads):

		this_read = bam_reads[read_upto]
		this_read_count = bam_reads_count[read_upto]
		this_read_compare = bam_reads_compare[read_upto]
		this_read_compare_collapsed = collapse_seq(this_read_compare)

		continue_looking_at_positions = True
		pos_upto = sliding_window - 1
		if (pos_upto >= len(this_read)):
			continue_looking_at_positions = False
		found_matching_region = False
		found_mismatching_region = False
		rightmost_matching_region = ''
		leftmost_mismatching_region =''
		# Start from the left, look for the right-most end of the matching region, having mismatching region to the right of it. That is the mobile element insert point.
		while (continue_looking_at_positions):

			this_read_pos_bp = this_read[pos_upto]
			this_ref_pos_bp = ref_seq[pos_upto]

			# if ( (found_matching_region == False) or ((found_matching_region == True) and (found_mismatching_region == False)) ):
			if (found_mismatching_region == False):
				if (this_read_pos_bp == this_ref_pos_bp):
					pos1 = pos_upto - sliding_window + 1
					pos2 = pos_upto + 1
					num_match = this_read_compare_collapsed.count( '1', pos1, pos2 ) # compare the sliding window to the left of and including up to pos_upto
					if (num_match >= minimum_match_amount):
						found_matching_region = True
						rightmost_matching_region = pos_upto

			if (found_matching_region):
				next_pos_upto = pos_upto + 1
				if (next_pos_upto < len(this_read)):
					next_read_pos_bp = this_read[next_pos_upto]
					next_ref_pos_bp = ref_seq[next_pos_upto]
					if ((next_read_pos_bp != next_ref_pos_bp) and (next_read_pos_bp != ' ') and (next_ref_pos_bp != ' ')):
						pos1 = pos_upto + 1
						pos2 = pos_upto + sliding_window
						if (pos2 < len(this_read)):
							num_mismatch = this_read_compare_collapsed.count( '0', pos1, pos2 ) # compare the sliding window to the right of pos_upto starting at the position after pos_upto			
							if (num_mismatch >= minimum_mismatch_amount):
								found_mismatching_region = True
								leftmost_mismatching_region = next_pos_upto

			if ((found_mismatching_region) and (found_matching_region)):

				this_seq = collapse_seq(this_read)
				this_seq_stripped = this_seq.strip()
				this_seq_start_pos = this_seq.find( this_seq_stripped )
				pos2 = pos_upto + 1
				pos1 = pos_upto - 50 # don't use pos1 = 0 because that makes this_percent_match_on_the_right dependent on concatenated read length, not on how well it matches to the left of this position 

				this_average_counts_on_the_right = -1
				this_mismatch_region_length_on_right = -1
				this_percent_match_on_the_left = -1

				# If there is more than one 5-prime end in the vicinity, 
				# then go for the one whose border3 position is after the border5 position of the 3-prime end found previously.
				# pos_upto and the_3prime_end are zero-based index

				choose_this_sequence = False
				if (best_pos < the_3prime_end):
					if (pos_upto > the_3prime_end): # (pos_upto > the_3prime_end) and (best_pos < the_3prime_end) so this pos_upto is clearly better than best_pos
						choose_this_sequence = True
					else: # pos_upto is also < the_3prime_end, as is best_pos. We would prefer a best_pos that is after the_3prime_end
						this_diff = abs(the_3prime_end - pos_upto)
						best_diff = abs(the_3prime_end - best_pos)
						if (this_diff < best_diff): # pick the border3 5prime-end that is closest to border5 3prime-end, unless the distances are almost the same (within 10% of each other)
							percent_increase_in_distance_from_previously_found_border = 0
							if (best_diff > 0):
								percent_increase_in_distance_from_previously_found_border = abs(this_diff - best_diff) * 100 / best_diff
							if (percent_increase_in_distance_from_previously_found_border > 10):
								choose_this_sequence = True
							else: # (percent_increase_in_distance_from_previously_found_border < 10)
								this_average_counts_on_the_right = calculate_average_counts( this_read_count[ (pos_upto+1) : len(this_read) ] )
								# if the difference in distance is only 10% or less, then choose based on number of supporting reads or else length of mismatch region
								if ((this_average_counts_on_the_right > 1) and (best_average_counts_on_the_right <= 1)):
									choose_this_sequence = True
								elif ((this_average_counts_on_the_right <= 1) and (best_average_counts_on_the_right > 1)):
									choose_this_sequence = False
								else: # they both have > 1 supporting read, or they both have only 1 supporting read, so choose the longest mismatch region
									this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )
									if (this_mismatch_region_length_on_right > best_mismatch_region_length_on_right):
										choose_this_sequence = True
				else: # best_pos > the_3prime_end. We want a best_pos > the_3prime_end.
					if (pos_upto > the_3prime_end): # both (best_pos > the_3prime_end) and (pos_upto > the_3prime_end), both are good
						this_diff = abs(pos_upto - the_3prime_end)
						best_diff = abs(best_pos - the_3prime_end)
						if ((this_diff <= 100) and (best_diff <= 100)): # both are close to the previously found border, so both are very good, choose based on mismatch_region, then supporting read counts, then matching ref.seq.

							# how long is this mismatching region on the right of this position
							this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )

							if ((this_mismatch_region_length_on_right >= 10) and (best_mismatch_region_length_on_right < 10)):
								choose_this_sequence = True
							elif ((this_mismatch_region_length_on_right < 10) and (best_mismatch_region_length_on_right >= 10)):
								choose_this_sequence = False
							elif ((this_mismatch_region_length_on_right >= 10) and (best_mismatch_region_length_on_right >= 10)): 
								# both sequences are good representations of this MEI when considering only length of mismatch_region_length_on_right

								this_average_counts_on_the_right = calculate_average_counts( this_read_count[ (pos_upto+1) : len(this_read) ] )
								# if the difference in distance is only 10% or less, then choose based on number of supporting reads or else length of mismatch region
								if ((this_average_counts_on_the_right > 1) and (best_average_counts_on_the_right <= 1)):
									choose_this_sequence = True
								elif ((this_average_counts_on_the_right <= 1) and (best_average_counts_on_the_right > 1)):
									choose_this_sequence = False
								else: # they both have > 1 supporting read, or they both have only 1 supporting read, so choose via another criteria

									# how much increase in this mismatching region on the right of this position
									percent_increase_in_mismatch_region_on_right = 0
									if (best_pos == -1):
										percent_increase_in_mismatch_region_on_right = 100
									else:
										if ((best_mismatch_region_length_on_right == 0) and (this_mismatch_region_length_on_right > 0)):
											percent_increase_in_mismatch_region_on_right = 100
										else:
											if (this_mismatch_region_length_on_right > best_mismatch_region_length_on_right):
												percent_increase_in_mismatch_region_on_right = 0
												if (best_mismatch_region_length_on_right > 0):
													percent_increase_in_mismatch_region_on_right = 0
													if (best_mismatch_region_length_on_right > 0):
														percent_increase_in_mismatch_region_on_right = ((this_mismatch_region_length_on_right - best_mismatch_region_length_on_right) * 100) / best_mismatch_region_length_on_right

									if (percent_increase_in_mismatch_region_on_right >= 10):
										choose_this_sequence = True
									else: # (percent_increase_in_mismatch_region_on_right < 10%) so we consider the 2 values to be roughly equal

										if (this_average_counts_on_the_right > best_average_counts_on_the_right):
											choose_this_sequence = True
										elif (this_mismatch_region_length_on_right > best_mismatch_region_length_on_right):
											choose_this_sequence = True
										else:

											# how well does this read match the ref.seq. on the left before this pos_upto after which it doesn't match
											this_percent_match_on_the_left = 0
											denom = pos2 - pos1 + 1
											if (denom > 0):
												this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom

											if (this_percent_match_on_the_left > best_percent_match_on_the_left):
												choose_this_sequence = True
											elif (best_average_counts_on_the_right > this_average_counts_on_the_right):
												choose_this_sequence = False
											elif (best_mismatch_region_length_on_right > this_mismatch_region_length_on_right):
												choose_this_sequence = False
											elif (best_percent_match_on_the_left > this_percent_match_on_the_left):
												choose_this_sequence = False
											else: # the 2 sequences are just as good or bad as each other as representing this MEI
												choose_this_sequence = False

							else: # ((this_mismatch_region_length_on_right < 10 base-pairs) and (best_mismatch_region_length_on_right < 10 base-pairs))
								# both sequences are poor representations of this MEI when considering only length of mismatch_region_length_on_right
								# however, they both have more than one supporting read

								if (this_average_counts_on_the_right > best_average_counts_on_the_right):
									choose_this_sequence = True
								elif (this_mismatch_region_length_on_right > best_mismatch_region_length_on_right):
									choose_this_sequence = True
								else:

									# how well does this read match the ref.seq. on the left before this pos_upto after which it doesn't match
									this_percent_match_on_the_left = 0
									denom = pos2 - pos1 + 1
									if (denom > 0):
										this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom

									if (this_percent_match_on_the_left > best_percent_match_on_the_left):
										choose_this_sequence = True
									elif (best_average_counts_on_the_right > this_average_counts_on_the_right):
										choose_this_sequence = False
									elif (best_mismatch_region_length_on_right > this_mismatch_region_length_on_right):
										choose_this_sequence = False
									elif (best_percent_match_on_the_left > this_percent_match_on_the_left):
										choose_this_sequence = False
									else: # the 2 sequences are just as good or bad as each other as representing this MEI
										choose_this_sequence = False

						else:
							if (this_diff < best_diff): # pick the border3 5prime-end that is closest to border5 3prime-end, unless the distances are almost the same (within 10% of each other)
								percent_increase_in_distance_from_previously_found_border = 0
								if (best_diff > 0):
									percent_increase_in_distance_from_previously_found_border = abs(this_diff - best_diff) * 100 / best_diff
								if (percent_increase_in_distance_from_previously_found_border > 10):
									choose_this_sequence = True
								else: # (percent_increase_in_distance_from_previously_found_border < 10)

									this_average_counts_on_the_right = calculate_average_counts( this_read_count[ (pos_upto+1) : len(this_read) ] )
									# how much increase in the number of reads showing this alternate sequence on the left
									percent_increase_in_average_counts_on_the_right = 0
									if (best_pos == -1):
										percent_increase_in_average_counts_on_the_right = 100
									else:
										if ((best_average_counts_on_the_right == 0) and (this_average_counts_on_the_right > 0)):
											percent_increase_in_average_counts_on_the_right = 100
										else:
											if (this_average_counts_on_the_right > best_average_counts_on_the_right):
												percent_increase_in_average_counts_on_the_right = 0
												if (best_average_counts_on_the_right > 0):
													percent_increase_in_average_counts_on_the_right = ((this_average_counts_on_the_right - best_average_counts_on_the_right) * 100) / best_average_counts_on_the_right

									if (percent_increase_in_average_counts_on_the_right >= 10):
										choose_this_sequence = True
									else: # (percent_increase_in_average_counts_on_the_right < 10%) so we consider the 2 values to be roughly equal

										# if the difference in distance is only 10% or less, then choose based on number of supporting reads or else length of mismatch region
										if ((this_average_counts_on_the_right > 1) and (best_average_counts_on_the_right <= 1)):
											choose_this_sequence = True
										elif ((this_average_counts_on_the_right <= 1) and (best_average_counts_on_the_right > 1)):
											choose_this_sequence = False
										else: # they both have > 1 supporting read, or they both have only 1 supporting read, so choose the longest mismatch region
											this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )
											if (this_mismatch_region_length_on_right > best_mismatch_region_length_on_right):
												choose_this_sequence = True
					# else (pos_upto <= the_3prime_end) and (best_pos > the_3prime_end), so keep the current best_pos, it is better than this pos_upto

				if (choose_this_sequence):
					if (this_average_counts_on_the_right == -1):
						this_average_counts_on_the_right = calculate_average_counts( this_read_count[ (pos_upto+1) : len(this_read) ] )
					if (this_mismatch_region_length_on_right == -1):
						this_mismatch_region_length_on_right = calculate_mismatch_region_length_from_the_left_end( this_read_compare_collapsed[ (pos_upto+1) : len(this_read) ], minimum_mismatch_percent )
					if (this_percent_match_on_the_left == -1):
						this_percent_match_on_the_left = 0
						denom = pos2 - pos1 + 1
						if (denom > 0):
							this_percent_match_on_the_left = this_read_compare_collapsed.count( '1', pos1, pos2 ) * 100 / denom
					best_pos = pos_upto
					best_seq = this_seq
					best_seq_stripped = this_seq_stripped
					best_seq_start_pos = this_seq_start_pos
					best_percent_match_on_the_left = this_percent_match_on_the_left
					best_average_counts_on_the_right = this_average_counts_on_the_right
					best_mismatch_region_length_on_right = this_mismatch_region_length_on_right

				# We found a possible MEI position in this sequence, whether or not we decided to use it.
				# Continue looking along this sequence in case there are even better possible MEI positions in this same sequence.
				# (If we don't do this, then when this sequence contains a strong MEI further along, 
				# the current weak MEI will be replaced by some other sequence's weak MEI, and we'll miss the strong MEI.)
				pos_upto = pos_upto + 1
				found_matching_region = False
				found_mismatching_region = False
				rightmost_matching_region = ''
				leftmost_mismatching_region =''
				if (pos_upto >= len(this_read)):
					continue_looking_at_positions = False

			else:
				pos_upto = pos_upto + 1
				if (pos_upto >= len(this_read)):
					continue_looking_at_positions = False

		read_upto = read_upto + 1
		if (read_upto >= len(bam_reads)):
			continue_looking_at_reads = False

	if (best_pos != -1):
		return_result = str(best_pos) + ',' + str(best_seq_start_pos) + ',' + str(best_average_counts_on_the_right) + ',' + str(best_mismatch_region_length_on_right) + ',' + str(best_percent_match_on_the_left) + ',' + best_seq_stripped
        return return_result

######################################################
def find_ref_seq_spanning_MEI_region_DO_BAM_READS_COVER_THIS_REGION( chrom, pos1_in_genome, pos2_in_genome, where_is_bam ):

	global ref_seq_startpos
	global ref_seq
	global read_length

	ref_seq_is_seen = False
	num_times_ref_seq_is_seen_in_this_region = 0
	pos1_in_ref_seq = pos1_in_genome - ref_seq_startpos # ref_seq_startpos is 1-based genome coordinate position and pos1_in_ref_seq is zero-based index
	# pos2_in_ref_seq = pos2_in_genome - ref_seq_startpos
	compare_length = pos2_in_genome - pos1_in_genome
	if (compare_length > 0):

		outlines = []
		samtools_command = 'samtools view ' + where_is_bam + ' ' + str(chrom) + ':' + str(pos1_in_genome) + '-' + str(pos2_in_genome)
		command_status, command_output = commands.getstatusoutput( samtools_command )
		if (command_status != 0):
			raise ValueError("\n\nIn Mobelwrap_convert_Mobster_predictions_into_Refined_predictions.py, was not able to read from BAM file using command:\n" + samtools_command + "\nThus will not continue processing any more sequences.\n")
		else:
			outlines = parse_out_any_warnings( command_output, samtools_command )

		for bam_record in outlines:
			bam_record = bam_record.strip()
			if (bam_record != ''):
				fields = bam_record.split("\t")
				# read_name = fields[0]
				# read_flag = fields[1]
				read_chrom = fields[2]
				read_pos = int(fields[3])
				# read_mapping_quality = int(fields[4])
				read_cigar = fields[5]
				#read_next_chrom = fields[6]
				# if (read_next_chrom == '='):
				#	read_next_chrom = read_chrom
				# read_next_pos = int(fields[7])
				# read_template_length = int(fields[8])
				read_sequence = fields[9]
				# read_sequence_quality = fields[10]
				read_end_pos = read_pos + len(read_sequence) - 1

				# compare this bam read to the ref.seq. for only the region we're interested in
				pos1_in_bam_read = pos1_in_genome - read_pos # zero-based index
				pos2_in_bam_read = pos2_in_genome - read_pos # this value is 1 past the end-pos of the zero-based index
				pos1_in_ref_seq = pos1_in_genome - ref_seq_startpos
				pos2_in_ref_seq = pos2_in_genome - ref_seq_startpos

				# The region of interest may partly fall outside the ref.seq. if we have found a border5 or border3 close to (within 20 bp) the edge of ref.seq. region. Adjust if necessary.
				if (pos1_in_ref_seq < 0):
					diff = 0 - pos1_in_ref_seq
					pos1_in_ref_seq = 0 # = pos1_in_ref_seq + diff
					pos1_in_bam_read = pos1_in_bam_read + diff
				if (pos2_in_ref_seq >= len(ref_seq)):
					diff = pos2_in_ref_seq - len(ref_seq) + 1
					pos2_in_ref_seq = pos2_in_ref_seq - diff
					pos2_in_bam_read = pos2_in_bam_read - diff

				# Only look at this read to see if it is the ref. seq. only if it can possibly contain the entire region of interest
				if ((pos1_in_bam_read >= 0) and (pos2_in_bam_read < len(read_sequence))):

					j = pos1_in_ref_seq
					count_matches = 0
					for i in range( pos1_in_bam_read, pos2_in_bam_read ):
						if ( read_sequence[i] == ref_seq[j] ): # if ( read_sequence[i:(i+1)] == ref_seq[j:(j+1)] ):
							count_matches = count_matches + 1
						j = j + 1
					percent_matches = 0
					percent_matches = count_matches * 100 / compare_length
					if (percent_matches > 90):
						ref_seq_is_seen = True
						num_times_ref_seq_is_seen_in_this_region = num_times_ref_seq_is_seen_in_this_region + 1

					else:
						# is it a pattern such as 14M1D136M
						is_a_stutter_overlapping_our_region = False
						regex1 = re.compile(r'[0-9]+M[0-9]D[0-9]+M$')
						regex_result = re.match( regex1, read_cigar )
						if (regex_result):
							split1 = read_cigar.split("M")
							num1 = int(split1[0])
							bit1 = split1[1]
							split2 = bit1.split("D")
							num2 = int(split2[0])
							num3 = int(split2[1])
							if ((num1 > 5) and (num2 <= 5) and (num3 > 5)):
								if ( ((read_pos <= pos1_in_genome) and (read_end_pos >= pos1_in_genome)) or ((read_pos <= pos2_in_genome) and (read_end_pos >= pos2_in_genome)) ):
									is_a_stutter_overlapping_our_region = True

						# is it a pattern such as 100M1I48M
						if (is_a_stutter_overlapping_our_region == False):
							regex2 = re.compile(r'[0-9]+M[0-9]I[0-9]+M$')
							regex_result = re.match( regex2, read_cigar )
							if (regex_result):
								split1 = read_cigar.split("M")
								num1 = int(split1[0])
								bit1 = split1[1]
								split2 = bit1.split("I")
								num2 = int(split2[0])
								num3 = int(split2[1])
								if ((num1 > 5) and (num2 <= 5) and (num3 > 5)):
									if ( ((read_pos <= pos1_in_genome) and (read_end_pos >= pos1_in_genome)) or ((read_pos <= pos2_in_genome) and (read_end_pos >= pos2_in_genome)) ):
										is_a_stutter_overlapping_our_region = True

						# This bam read is a stutter match of ref.seq. so let's count it
						if (is_a_stutter_overlapping_our_region):
							ref_seq_is_seen = True
							num_times_ref_seq_is_seen_in_this_region = num_times_ref_seq_is_seen_in_this_region + 1

        return ref_seq_is_seen, num_times_ref_seq_is_seen_in_this_region

######################################################
def find_ref_seq_spanning_MEI_region( chrom, MEI_pos1_in_genome, MEI_pos2_in_genome, where_is_bam ):

	# Read the bam file to see if there are any reads spanning this region of a new MEI, that match the ref.seq.
	# We read the bamfile again even though we already have concatenated bam reads from identify_split_read_clusters.py
	# because identify_split_read_clusters.py concatenates when there is overlapping sequence
	# and does concatenate sequences from both sides of an MEI even when there is no read containing both sides of the MEI.

	# Look for ref.seq. spanning 50 base-pairs on either side of border5, and 50 base-pairs on either side of border3.
	# If we see either one, then we have seen the ref.seq. on at least one of the borders and thus it is not a homozygous MEI.

	# Heterozygous will return 'Seen' to say that ref.seq. was seen, and will become genotype 0/1 downstream.
	# Homozygous will return 'NotSeen' to say that ref.seq. was not seen, and will become genotype 1/1 downstream.
	# However, if border5 and border3 are further apart than the read_length, then we can't tell whether it is homozygous or heterozygous
	# because the reads are not long enough to have both ends of the MEI on the one read.
	# 'SeenUndetermined' will become genotype 1/. downstream, but we are always chosing Seen or NotSeen, so SeenUndetermined won't happen.

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq_startpos
	global ref_seq
	global read_length

	return_result = ''
	ref_seq_spanning_sequence_seen = 'SeenUndetermined'
	ref_seq_spanning_sequence_seen_total_reads_seen = 0

	if (MEI_pos1_in_genome > MEI_pos2_in_genome):
		MEI_pos2_in_genome_saved = MEI_pos2_in_genome
		MEI_pos2_in_genome = MEI_pos1_in_genome
		MEI_pos1_in_genome = MEI_pos2_in_genome_saved

	# Look for ref.seq. in region spanning 50 base-pairs either side of border5

	pos1_in_genome = MEI_pos1_in_genome - 50
	pos2_in_genome = MEI_pos1_in_genome + 50
	border5_ref_seq_seen = False
	border5_total_reads_seen = 0
	if ((pos2_in_genome - pos1_in_genome) <= read_length):

		ref_seq_is_seen, num_times_ref_seq_is_seen_in_this_region = find_ref_seq_spanning_MEI_region_DO_BAM_READS_COVER_THIS_REGION( chrom, pos1_in_genome, pos2_in_genome, where_is_bam )

		if (ref_seq_is_seen):
			border5_ref_seq_seen = 'Seen'
			border5_total_reads_seen = num_times_ref_seq_is_seen_in_this_region

	# Look for ref.seq. in region spanning 50 base-pairs either side of border3

	pos1_in_genome = MEI_pos2_in_genome - 50
	pos2_in_genome = MEI_pos2_in_genome + 50
	border3_ref_seq_seen = False
	border3_total_reads_seen = 0
	if ((pos2_in_genome - pos1_in_genome) <= read_length):

		ref_seq_is_seen, num_times_ref_seq_is_seen_in_this_region = find_ref_seq_spanning_MEI_region_DO_BAM_READS_COVER_THIS_REGION( chrom, pos1_in_genome, pos2_in_genome, where_is_bam )

		if (ref_seq_is_seen):
			border3_ref_seq_seen = 'Seen'
			border3_total_reads_seen = num_times_ref_seq_is_seen_in_this_region

	if (border5_ref_seq_seen or border3_ref_seq_seen):
		ref_seq_spanning_sequence_seen = 'Seen'
		ref_seq_spanning_sequence_seen_total_reads_seen = border5_total_reads_seen + border3_total_reads_seen
	else:
		ref_seq_spanning_sequence_seen = 'NotSeen'

        return ref_seq_spanning_sequence_seen, ref_seq_spanning_sequence_seen_total_reads_seen

######################################################
def is_this_ME_in_existing_ME_region( chrom, pos1, pos2, mobster_repmask_file_name ):

	global existing_ME_regions_chrom
	global existing_ME_regions_border5
	global existing_ME_regions_border3
	global args

	# if we haven't already read in this information, then do so
	if (len(existing_ME_regions_border5) == 0):
		input_repmask = open(mobster_repmask_file_name, 'r')
		in_header = True
		for inline in input_repmask:
			inline = inline.strip()
			if (inline != ''):
				if (in_header):
					in_header = False
				else:
					infields = inline.split("\t")
					genoName = infields[5]
					remove_chr = False
					if args.keep_chr_prefix is not None:
						remove_chr = True
						if (args.keep_chr_prefix == 'SOME'):
							if (len(genoName) >= 4):
								chr4 = genoName[3:4]
								if (chr4 in ('1','2','3','4','5','6','7','8','9','X','Y','M')):
									remove_chr = True
								else:
									remove_chr = False
					if (remove_chr):
						if (len(genoName) >= 3):
							if (genoName[0:3] == 'chr'):
								genoName = genoName[3:]
					genoName = str(genoName)
					genoStart = int(infields[6])
					genoEnd = int(infields[7])
					existing_ME_regions_chrom.append( genoName )
					existing_ME_regions_border5.append( genoStart )
					existing_ME_regions_border3.append( genoEnd )

	# is this ME in an existing region? 
	# Is one of it's borders within 90 base-pair of the border of an existing mobile element?
	is_in_existing_ME_region = False
	look_for_ME = True
	i = 0
	while (look_for_ME):
		if (chrom == existing_ME_regions_chrom[i]):
			diff1 = abs(pos1 - existing_ME_regions_border5[i])
			diff2 = abs(pos1 - existing_ME_regions_border3[i])
			diff3 = abs(pos2 - existing_ME_regions_border5[i])
			diff4 = abs(pos2 - existing_ME_regions_border3[i])
			diff_min = min( diff1, diff2, diff3, diff4 )
			if (diff_min <= 90):
				is_in_existing_ME_region = True
				look_for_ME = False
		i = i + 1
		if (i >= len(existing_ME_regions_border5)):
			look_for_ME = False

        return is_in_existing_ME_region

######################################################
def process_one_input_line_from_processing_batch( processing_batch_input_line ):

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq
	global ref_seq_startpos
	global read_length
	global existing_ME_regions_chrom
	global existing_ME_regions_border5
	global existing_ME_regions_border3
	global args
	global python_path_2

	return_array = [ '', '' ]

	read_length = int(args.read_length_of_in_bam_reads)
	str_read_length = str(read_length)

	inline = processing_batch_input_line
	if (inline != ''):

		wrote_out_this_record = False
		inline = inline.replace( "\n", '' )
		inline = inline.replace( "\r", '' )
		if (inline != ''):

			infields = inline.split("\t")
			chrom = infields[0] # Chr
			remove_chr = False
			if args.keep_chr_prefix is not None:
				remove_chr = True
				if (args.keep_chr_prefix == 'SOME'):
					if (len(chrom) >= 4):
						chr4 = chrom[3:4]
						if (chr4 in ('1','2','3','4','5','6','7','8','9','X','Y','M')):
							remove_chr = True
						else:
							remove_chr = False
			if (remove_chr):
				# Mobster puts a 'chr' in front of chromosome. Remove it
				if (len(chrom) >= 3):
					if (chrom[0:3] == 'chr'):
						chrom = chrom[3:]
						inline = inline[3:]
			chrom = str(chrom)

			# ignore reads mapped to the following chromosomes:
			# NC_007605	:	herpes virus
			# hs37d5	:	decoy
			# NC_001422	:	phiX, used as an Illumina control
			# *		:	to be ignored

			if ((chrom[0:9] == 'NC_007605') or (chrom == 'hs37d5') or (chrom[0:9] == 'NC_001422') or (chrom == '*')):

				ignore_this_sequence = True
				if args.discard_Mobster_file is not None:
					qual = ''
					outline = 'Ignore_this_contig_or_chrom' + "\t" + inline + "\t" + qual + "\n"
					return_array[1] = outline

			else:

				mobile_element_type = str(infields[1])
				insert_point = int(infields[2])
				border5 = int(infields[3])
				border3 = int(infields[4])
				# merged_field = int(infields[5])
				# sample_field = int(infields[6])
				# sample_counts = int(infields[7])
				# cluster5_length = convert_to_integer(infields[8])
				# cluster3_length = convert_to_integer(infields[9])
				cluster5_hits = convert_to_integer(infields[10])
				cluster3_hits = convert_to_integer(infields[11])
				split5_hits = convert_to_integer(infields[12])
				split3_hits = convert_to_integer(infields[13])
				polyA5_hits = convert_to_integer(infields[14])
				polyT5_hits = convert_to_integer(infields[15])
				polyA3_hits = convert_to_integer(infields[16])
				polyT3_hits = convert_to_integer(infields[17])
				original_discordant_unique = convert_to_integer( infields[18] )
				if (insert_point < 1):
					insert_point = 1
				if (border5 < 1):
					border5 = 1
				if (border3 < 1):
					border3 = 1
				pos1 = min( insert_point, border5, border3 )
				pos2 = max( insert_point, border5, border3 )
				pos1 = str(pos1)
				pos2 = str(pos2)

				# Calculate the quality of this MEI call from Mobster counts.
				# If necessary, filter this MEI call by quality
				discord5 = cluster5_hits - split5_hits
				discord3 = cluster3_hits - split3_hits
				Qsplit5 = for_qual_calc( split5_hits )
				Qsplit3 = for_qual_calc( split3_hits )
				Qdiscord5 = for_qual_calc( discord5 )
				Qdiscord3 = for_qual_calc( discord3 )
				QpolyA5 = for_qual_calc( polyA5_hits )
				QpolyT5 = for_qual_calc( polyT5_hits )
				QpolyA3 = for_qual_calc( polyA3_hits )
				QpolyT3 = for_qual_calc( polyT3_hits )
				Qoriginal_discordant_unique = for_qual_calc( original_discordant_unique )
				qual = math.log( (Qsplit5*10 + Qsplit3*10 + Qdiscord5*2 + Qdiscord3*2 + QpolyA5*2 + QpolyT5*2 + QpolyA3*2 + QpolyT3*2 + Qoriginal_discordant_unique*2), 10 ) * 10
				qual = str(int(round(qual)))

				discard_this_MEI_call = False
				if args.quality_filter is not None:
					if (int(qual) < int(args.quality_filter)):
						discard_this_MEI_call = True

				if (discard_this_MEI_call):

					ignore_this_sequence = True
					if args.discard_Mobster_file is not None:
						outline = 'Quality_score_below_filter_level' + "\t" + inline + "\t" + qual + "\n"
						return_array[1] = outline

				else:

					# if ((split5_hits > 0) or (split3_hits > 0) or (polyA5_hits > 0) or (polyT5_hits > 0) or (polyA3_hits > 0) or (polyT3_hits > 0)):

					# get the list of unique concatenated BAM reads for this mobile-element-insertion's region

					new_bam_reads = []
					new_bam_reads_count = []
					bam_reads_compare = []
					bam_cluster_position = chrom + ':' + pos1 + '-' + pos2
					bam_cluster_command = 'python ' + python_path_2 + 'Mobelwrap_identify_split_read_clusters.py ' + bam_cluster_position + ' ' + args.in_bam_reads + ' ' + str_read_length + ' ' + args.in_reference_genome
					command_status, command_output = commands.getstatusoutput( bam_cluster_command )
					if (command_status != 0):
						raise ValueError("\n\nWas not able to read the input BAM file via python program " + bam_cluster_command + "\nReceived status: " + str(command_status) + "\nThus will not continue processing any more mobile elements.\n\n")
					command_output_lines = command_output.split("\n")
					bam_reads_chrom = command_output_lines[0]
					ref_seq_startpos = int(command_output_lines[1])
					for i in range( 2, len(command_output_lines), 2 ):
						if (command_output_lines[i] != ''):
							bam_read_seq = command_output_lines[i]
							bam_read_seq_count = command_output_lines[i+1]
							bam_read = []
							bam_read_count = []
							bam_read = [' '] * len(bam_read_seq)
							bam_read_count = [0] * len(bam_read_seq)
							for j in range( 0, len(bam_read_seq) ):
								this_bp = bam_read_seq[j]
								bam_read[j] = this_bp
								this_bp_count = bam_read_seq_count[j]
								bam_read_count[j] = this_bp_count
							new_bam_reads.append( bam_read )
							new_bam_reads_count.append( bam_read_count )
					bam_reads = new_bam_reads
					bam_reads_count = new_bam_reads_count
					bam_read_length = 0
					new_pos1 = pos1
					new_pos2 = pos1
					if (len(bam_reads) > 0):

						# Can continue looking at the concatenated bam reads only if there are some.

						first_bam_read = bam_reads[0]
						bam_read_length = len(first_bam_read)
						new_pos1 = ref_seq_startpos
						new_pos2 = new_pos1 + bam_read_length - 1

						# get the reference sequence for the region of this mobile-element-insertion

						command_output_lines = []
						ref_seq = ''
						samtools_faidx_position = chrom + ':' + str(new_pos1) + '-' + str(new_pos2)
						samtools_faidx_command = 'samtools faidx ' + args.in_reference_genome + ' ' + samtools_faidx_position
						command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
						if (command_status != 0):
							raise ValueError("\n\nIn Mobelwrap_convert_Mobster_predictions_into_Refined_predictions.py, was not able to get the reference sequence from reference genome for this mobile element's region using command:\n" + samtools_faidx_command + "\nReceived status: " + str(command_status) + "\nThus will not continue processing any more mobile elements.\n")
						else:
							command_output_lines = parse_out_any_warnings( command_output, samtools_faidx_command )

						# some positions may not be found in reference genome, will have a header but no sequence
						if (len(command_output_lines) > 1):
							ref_seq_string = ''
							for i in range( 1, len(command_output_lines) ):
								ref_seq_string = ref_seq_string + command_output_lines[i]
							for i in range( 0, len(ref_seq_string) ):
								this_bp = ref_seq_string[i:i+1]
								ref_seq = ref_seq + this_bp
						ref_seq_expected_length = new_pos2 - new_pos1 + 1
						ref_seq_actual_length = len(ref_seq)
						if (ref_seq_actual_length < ref_seq_expected_length):
							ref_seq = ref_seq + (' ' * (ref_seq_expected_length - ref_seq_actual_length))

						# For each concatenated bam read, figure out how similar it is to the ref.seq.

						continue_looking_at_reads = True
						read_upto = 0
						if (read_upto >= len(bam_reads)):
							continue_looking_at_reads = False
						while (continue_looking_at_reads):
							this_read = bam_reads[read_upto]
							this_read_compare = [' '] * len(ref_seq)
							smallest_read_length = min( len(ref_seq), len(this_read) )
							for i in range( 0, smallest_read_length ):
								if ((ref_seq[i] != ' ') and (ref_seq[i] != 'N') and (this_read[i] != ' ') and (this_read[i] != 'N')):
									if (ref_seq[i] == this_read[i]):
										this_read_compare[i] = '1'
									else:
										this_read_compare[i] = '0'
							bam_reads_compare.append( this_read_compare )
							read_upto = read_upto + 1
							if (read_upto >= len(bam_reads)):
								continue_looking_at_reads = False

						# Find the 3-prime end and 5-prime end of this mobile-element-insertion (MEI) in the unique concatenated BAM reads.
						# The position returned is the last position matching the reference genome before the unmatching mobile element sequence.
						# In these comments, 3-prime refers to the 3-prime end of the MEI (not the 3-prime of the reference sequence).
						# In these comments, 5-prime refers to the 5-prime end of the MEI (not the 5-prime of the reference sequence).
						# Find the 3-prime end first. 
						# If there are multiple possible 3-prime end sequences, 
						# choose the one with the most amounts of reads in the part of the read that does not match the reference sequence.
						# That will be the most obvious-looking insertion.
						# If there is more than one 5-prime end in the vicinity, then go for the one that is on the same read as the 3-prime end.
						# If there isn't one, then go for the one that whose position is the closest to the position of the 3-prime end.
						# Otherwise, when finding a 5-prime end, favour a 5-prime end whose position is after the position of the 3-prime end.
						# Why? Because that is the pattern for mobile element insertions that we saw in MGRB data whilst developing this program.
						# There was an example where the 5-prime position after the 3-prime position was the same mobile element insertion,
						# whilst the 5-prime position before the 3-prime position in fact seemed to belong to the preceding mobile element insertion.

						the_3prime_end_pos = -1
						the_5prime_end_pos = -1
						the_3prime_end_seq_start_pos = -1
						the_5prime_end_seq_start_pos = -1
						the_3prime_end_average_read_count_for_MEI_read_fragment = -1
						the_5prime_end_average_read_count_for_MEI_read_fragment = -1
						the_3prime_end_length_of_MEI_read_fragment = -1
						the_5prime_end_length_of_MEI_read_fragment = -1
						the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = -1
						the_5prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = -1
						the_3prime_end_seq = ''
						the_5prime_end_seq = ''
						inhouse_ref_seq_spanning_sequence_seen = ''
						inhouse_ref_seq_spanning_sequence_seen_start_pos = ''
						inhouse_ref_seq_spanning_sequence_seen_end_pos = ''
						inhouse_ref_seq_spanning_sequence_seen_num_reads = ''

						the_3prime_end = find_3prime_end_of_mobile_element_insertion() ### <=== look for border5, for an MEI insertion evidenced by a split read

						feed_3prime_end_to_5prime_search = -1
						if (the_3prime_end != ''):
							the_3prime_end_fields = the_3prime_end.split(',')
							the_3prime_end_pos = int(the_3prime_end_fields[0]) + ref_seq_startpos # the_3prime_end_fields[0] is zero-based index, so just add to ref_seq_startpos to get 1-based genome coordinate position
							the_3prime_end_seq_start_pos = int(the_3prime_end_fields[1]) + ref_seq_startpos
							the_3prime_end_average_read_count_for_MEI_read_fragment = float(the_3prime_end_fields[2])
							the_3prime_end_length_of_MEI_read_fragment = int(the_3prime_end_fields[3])
							the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = int(the_3prime_end_fields[4])
							the_3prime_end_seq = the_3prime_end_fields[5]
							feed_3prime_end_to_5prime_search = int(the_3prime_end_fields[0])

							the_5prime_end = find_5prime_end_of_mobile_element_insertion_got_3prime_already( feed_3prime_end_to_5prime_search ) ### <=== look for border3, for an MEI insertion evidenced by a split read

							if (the_5prime_end != ''):
								the_5prime_end_fields = the_5prime_end.split(',')
								the_5prime_end_pos = int(the_5prime_end_fields[0]) + ref_seq_startpos # the_5prime_end_fields[0] is zero-based index, so just add to ref_seq_startpos to get 1-based genome coordinate position
								the_5prime_end_seq_start_pos = int(the_5prime_end_fields[1]) + ref_seq_startpos
								the_5prime_end_average_read_count_for_MEI_read_fragment = float(the_5prime_end_fields[2])
								the_5prime_end_length_of_MEI_read_fragment = int(the_5prime_end_fields[3])
								the_5prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = int(the_5prime_end_fields[4])
								the_5prime_end_seq = the_5prime_end_fields[5]

								# If border3 is good and border5 is terrible, have another look for a better border5

								if ((the_3prime_end_pos != -1) and (the_5prime_end_pos != -1)):
									if (the_3prime_end_pos > the_5prime_end_pos): # this is border5 > border3, this is not good, we are more convinced it is a true MEI if border5 < border3
										if (the_3prime_end_average_read_count_for_MEI_read_fragment < the_5prime_end_average_read_count_for_MEI_read_fragment): 

											feed_5prime_end_to_3prime_search = int(the_5prime_end_fields[0])

											the_3prime_end = find_3prime_end_of_mobile_element_insertion_got_5prime_already( feed_5prime_end_to_3prime_search ) ### <=== look for border5, for an MEI insertion evidenced by a split read

											if (the_3prime_end != ''):
												the_3prime_end_fields = the_3prime_end.split(',')
												the_3prime_end_pos = int(the_3prime_end_fields[0]) + ref_seq_startpos # the_3prime_end_fields[0] is zero-based index, so just add to ref_seq_startpos to get 1-based genome coordinate position
												the_3prime_end_seq_start_pos = int(the_3prime_end_fields[1]) + ref_seq_startpos
												the_3prime_end_average_read_count_for_MEI_read_fragment = float(the_3prime_end_fields[2])
												the_3prime_end_length_of_MEI_read_fragment = int(the_3prime_end_fields[3])
												the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = int(the_3prime_end_fields[4])
												the_3prime_end_seq = the_3prime_end_fields[5]

						else:
							the_5prime_end = find_5prime_end_of_mobile_element_insertion() ### <=== look for border3, for an MEI insertion evidenced by a split read

							feed_5prime_end_to_3prime_search = -1
							if (the_5prime_end != ''):
								the_5prime_end_fields = the_5prime_end.split(',')
								the_5prime_end_pos = int(the_5prime_end_fields[0]) + ref_seq_startpos # the_5prime_end_fields[0] is zero-based index, so just add to ref_seq_startpos to get 1-based genome coordinate position
								the_5prime_end_seq_start_pos = int(the_5prime_end_fields[1]) + ref_seq_startpos
								the_5prime_end_average_read_count_for_MEI_read_fragment = float(the_5prime_end_fields[2])
								the_5prime_end_length_of_MEI_read_fragment = int(the_5prime_end_fields[3])
								the_5prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = int(the_5prime_end_fields[4])
								the_5prime_end_seq = the_5prime_end_fields[5]
								feed_5prime_end_to_3prime_search = int(the_5prime_end_fields[0])

								the_3prime_end = find_3prime_end_of_mobile_element_insertion_got_5prime_already( feed_5prime_end_to_3prime_search ) ### <=== look for border5, for an MEI insertion evidenced by a split read

								if (the_3prime_end != ''):
									the_3prime_end_fields = the_3prime_end.split(',')
									the_3prime_end_pos = int(the_3prime_end_fields[0]) + ref_seq_startpos # the_3prime_end_fields[0] is zero-based index, so just add to ref_seq_startpos to get 1-based genome coordinate position
									the_3prime_end_seq_start_pos = int(the_3prime_end_fields[1]) + ref_seq_startpos
									the_3prime_end_average_read_count_for_MEI_read_fragment = float(the_3prime_end_fields[2])
									the_3prime_end_length_of_MEI_read_fragment = int(the_3prime_end_fields[3])
									the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = int(the_3prime_end_fields[4])
									the_3prime_end_seq = the_3prime_end_fields[5]

						if ((the_5prime_end != '') or (the_3prime_end != '')):

							if ((the_5prime_end != '') and (the_3prime_end != '')):

								# We have seen this MEI call. Is it on one strand or both strands.
								# Look to see if there is the reference sequence on either side of border5 or either side of border3 positions (and 50 base-pairs beyond each).
								# If we see either of these, then we assume that we have seen the ref.seq. here and that it is not a homozygous mobile-element insertion.

								inhouse_ref_seq_spanning_sequence_seen = '' # Values are Seen, NotSeen
								inhouse_ref_seq_spanning_sequence_seen_start_pos = the_5prime_end_pos
								inhouse_ref_seq_spanning_sequence_seen_end_pos = the_3prime_end_pos
								if (the_3prime_end < the_5prime_end):
									inhouse_ref_seq_spanning_sequence_seen_start_pos = the_3prime_end_pos
									inhouse_ref_seq_spanning_sequence_seen_end_pos = the_5prime_end_pos
								inhouse_ref_seq_spanning_sequence_seen_num_reads = ''

								inhouse_ref_seq_spanning_sequence_seen, inhouse_ref_seq_spanning_sequence_seen_num_reads = find_ref_seq_spanning_MEI_region( chrom, inhouse_ref_seq_spanning_sequence_seen_start_pos, inhouse_ref_seq_spanning_sequence_seen_end_pos, args.in_bam_reads )

							# Now prepare all the counts and statistics about this MEI. They will be written to output.

							if (the_5prime_end_pos == -1):
								the_5prime_end_pos = ''
							if (the_5prime_end_seq_start_pos == -1):
								the_5prime_end_seq_start_pos = ''
							if (the_5prime_end_average_read_count_for_MEI_read_fragment == -1):
								the_5prime_end_average_read_count_for_MEI_read_fragment = ''
							if (the_5prime_end_length_of_MEI_read_fragment == -1):
								the_5prime_end_length_of_MEI_read_fragment = ''
							if (the_5prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference == -1):
								the_5prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = ''
							if (the_3prime_end_pos == -1):
								the_3prime_end_pos = ''
							if (the_3prime_end_seq_start_pos == -1):
								the_3prime_end_seq_start_pos = ''
							if (the_3prime_end_average_read_count_for_MEI_read_fragment == -1):
								the_3prime_end_average_read_count_for_MEI_read_fragment = ''
							if (the_3prime_end_length_of_MEI_read_fragment == -1):
								the_3prime_end_length_of_MEI_read_fragment = ''
							if (the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference == -1):
								the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference = ''

							outline = inline + "\t" + str(qual)
							outline = outline + "\t" + str(the_3prime_end_pos) + "\t" + str(the_5prime_end_pos)
							outline = outline + "\t" + str(the_3prime_end_seq_start_pos) + "\t" + str(the_5prime_end_seq_start_pos)
							outline = outline + "\t" + str(the_3prime_end_average_read_count_for_MEI_read_fragment) + "\t" + str(the_5prime_end_average_read_count_for_MEI_read_fragment)
							outline = outline + "\t" + str(the_3prime_end_length_of_MEI_read_fragment) + "\t" + str(the_5prime_end_length_of_MEI_read_fragment)
							outline = outline + "\t" + str(the_3prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference) + "\t" + str(the_5prime_end_percent_that_non_MEI_part_of_MEI_reads_matches_reference)
							outline = outline + "\t" + inhouse_ref_seq_spanning_sequence_seen + "\t" + str(inhouse_ref_seq_spanning_sequence_seen_start_pos)
							outline = outline + "\t" + str(inhouse_ref_seq_spanning_sequence_seen_end_pos) + "\t" + str(inhouse_ref_seq_spanning_sequence_seen_num_reads)
							outline = outline + "\t" + the_3prime_end_seq + "\t" + the_5prime_end_seq + "\n"
							return_array[0] = outline = outline
							wrote_out_this_record = True

						# If Mobster is run with the repmask file hg19_alul1svaerv.txt that contains known mobile-element-insertions (MEI) in the human reference genome,
						# then any MEI in a sample that fall within or close to that region (within 90 bp of an annotated MEI of the same ME family) will be ignored and not reported by Mobster.
						# If Mobster is run with an empty repmask file, then there is a possibility that an MEI reported by Mobster could actually be an MEI existing in the reference genome and is not a novel MEI for this sample.
						# The design of this convert_imprecise_Mobster_predictions_into_precise.py program is that those reference-genome-MEI will not be removed but will appear with IMPRECISE border5 and IMPRECISE border3.
						# Mobster-called MEI for which there are no split-reads support and the support is discordant reads,
						# will also come out of this convert_imprecise_Mobster_predictions_into_precise.py program with IMPRECISE border5 and IMPRECISE border3.
						# Currently, this convert_imprecise_Mobster_predictions_into_precise.py program does not tell the difference between Mobster-called reference-genome-MEI and Mobster-called discordant-read-support-only-MEI.
						# Both are retained by this convert_imprecise_Mobster_predictions_into_precise.py program.

					if (wrote_out_this_record == False):
						# If we didn't write out this ME (because we couldn't find any BAM reads showing its border5 or border3), then write it out.
						# If it is in an ME region already existing in the reference genome, then don't write it out, drop it (the reason we couldn't find BAM reads supporting this novel ME is that there is not a novel ME here)
						if (is_this_ME_in_existing_ME_region( chrom, border5, border3, args.existing_mobile_element_regions )):

							ignore_this_sequence = True
							if args.discard_Mobster_file is not None:
								outline = 'In_existing_ME_region_and_weak_read_support' + "\t" + inline + "\t" + qual + "\n"
								return_array[1] = outline

						else: # (is_this_ME_in_existing_ME_region( chrom, border5, border3, args.existing_mobile_element_regions ) == False):
							outline = inline + "\n"
							return_array[0] = outline

	return return_array

######################################################
def main():

	global bam_reads
	global bam_reads_count
	global bam_reads_compare
	global ref_seq_startpos
	global ref_seq
	global read_length
	global existing_ME_regions_chrom
	global existing_ME_regions_border5
	global existing_ME_regions_border3
	global args
	global python_path_2

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed Mobster predictions file containing mobile element insertion (MEI) calls for a sample, read in the sample\'s BAM file and the genome reference data file. For each MEI call, read BAM to work out it\'s exact location. Output the MEI sequences fragments, revised positions, and whether we can precisely determine those positiions, along with the input Mobster data.')
	parser.add_argument('-i', action="store", dest="in_Mobster", required=True, help='Input file containing Mobster predictions')
	parser.add_argument('-b', action="store", dest="in_bam_reads", required=True, help='Input BAM file containing mapped reads. Needs to be indexed. BAM reads in the region of each Mobster mobile-element-insertion-prediction will be retrieved from this BAM file to determine the exact insertions points of the Mobster mobile-element-insertion-prediction')
	parser.add_argument('-l', action="store", dest="read_length_of_in_bam_reads", required=True, help='The general read length of the reads in the input BAM file')
	parser.add_argument('-r', action="store", dest="in_reference_genome", required=True, help='Reference genome. A fasta file indexed by samtools faidx. BAM reads will be compared to this reference to determine the exact insertion points for Mobster mobile-element-insertion-predictions')
	parser.add_argument('-e', action="store", dest="existing_mobile_element_regions", required=True, help='Mobster\'s repmask file \'hg19_alul1svaerv.txt\' consisting of 1 line header and 1 line per mobile element region already existing in the reference genome.')
	parser.add_argument('-o', action="store", dest="out_Mobster_file", required=True, help='Output file will contain refined Mobster mobile-element-insertion-predictions, and where possible, extra information appended to the record to precisely identify the insertion points. (And if this refinement program does not find extra information for Mobster MEIs in existing ME region, then the MEI is dropped.)')
	parser.add_argument('-s', action="store", dest="restart_Mobster_file", required=False, help='If this file is present, then keep the MEI calls from this file and continue processing input Mobster calls from the next input Mobster MEI after the last MEI call seen in this file.')
	parser.add_argument('-q', action="store", dest="quality_filter", required=False, help='If present, then filter MEI calls, keeping calls having this or higher quality score.')
	parser.add_argument('-d', action="store", dest="discard_Mobster_file", required=False, help='If this output file is present, then record all discarded Mobster MEI calls to this file.')
	parser.add_argument('-c', action="store", dest="num_cores", required=False, help='Number of concurrent processors for multiprocessor processing')
	parser.add_argument('-k', action="store", dest="keep_chr_prefix", required=False, help='Mobster puts chr in front of chromosome and this program removes it. If this flag is present then this program will not remove the chr prefix.')
	parser.add_argument('-p', action="store", dest="python_path_2", required=False, help='The path where the python program called by this program is found. Alternatively, set PYTHONPATH2 environmental variable to this value.')
	args = parser.parse_args()

	# read each mobile-element-insertion-prediction record in the input Mobster file, and in the restart file if it is present

	python_path_2 = ''
	if args.python_path_2 is not None:
		python_path_2 = args.python_path_2
	else:
		if ('PYTHONPATH2' in os.environ):
			python_path_2 = os.environ['PYTHONPATH2']
	if (python_path_2 != ''):
		python_path_2 = python_path_2 + '/'

	output_Mobster = open(args.out_Mobster_file, 'w')
	input_Mobster_inlines = []
	input_Mobster = open(args.in_Mobster, 'r')
	for inline in input_Mobster:
		input_Mobster_inlines.append( inline )

	restart_Mobster_inlines = []
	if args.restart_Mobster_file is not None:
		restart_Mobster = open(args.restart_Mobster_file, 'r')
		for inline in restart_Mobster:
			restart_Mobster_inlines.append( inline )

	if args.discard_Mobster_file is not None:
		discard_Mobster = open(args.discard_Mobster_file, 'w')

	# If there is a restart file, check that it matches the input file.
	if args.restart_Mobster_file is not None:
		if ( len(restart_Mobster_inlines) > len(input_Mobster_inlines) ):
			raise ValueError("\n\nThe restart Mobster file " + args.restart_Mobster_file + " (to be read on input and contains refined Mobster output calls from the last time this program was called) has more lines (" + str(len(restart_Mobster_inlines)) + " lines) than the input Mobster file " + args.in_Mobster + " (" + str(len(input_Mobster_inlines)) + " lines).\nThese must not be the correct input files and processing will not be carried out.\n")
		if ( len(restart_Mobster_inlines) < 4 ):
			raise ValueError("\n\nThe restart Mobster file " + args.restart_Mobster_file + " doesn't seem to have a full header (has only " + str(len(restart_Mobster_inlines)) + " lines) and so processing will not be carried out.\n")
		in_upto = 4
		for restart_upto in range( 4, len(restart_Mobster_inlines) ):
			restart_MEI_call = restart_Mobster_inlines[restart_upto]
			restart_fields = restart_MEI_call.split("\t")
			restart_0 = restart_fields[0]
			restart_1 = restart_fields[1]
			restart_2 = restart_fields[2]
			restart_3 = restart_fields[3]
			restart_4 = restart_fields[4]
			look_for_this_restart_line_in_the_input_file = True
			while (look_for_this_restart_line_in_the_input_file):
				in_MEI_call = input_Mobster_inlines[in_upto]
				if (in_MEI_call.strip() == ''):
					ignore_empty_line = True
				else:
					in_fields = in_MEI_call.split("\t")
					in_0 = in_fields[0]
					in_1 = in_fields[1]
					in_2 = in_fields[2]
					in_3 = in_fields[3]
					in_4 = in_fields[4]
					remove_chr = False
					if args.keep_chr_prefix is not None:
						remove_chr = True
						if (args.keep_chr_prefix == 'SOME'):
							if (len(in_0) >= 4):
								chr4 = in_0[3:4]
								if (chr4 in ('1','2','3','4','5','6','7','8','9','X','Y','M')):
									remove_chr = True
								else:
									remove_chr = False
					if (remove_chr):
						if (len(in_0) >= 3):
							if (in_0[0:3] == 'chr'):
								in_0 = in_0[3:]
					if ((in_0 == restart_0) and (in_1 == restart_1) and (in_2 == restart_2) and (in_3 == restart_3) and (in_4 == restart_4)):
						look_for_this_restart_line_in_the_input_file = False
				if (look_for_this_restart_line_in_the_input_file == True):
					in_upto = in_upto + 1
					if (in_upto >= len(input_Mobster_inlines)):
						line_upto = restart_upto + 1
						raise ValueError("\n\nLine " + str(line_upto) + " of the restart file " + args.restart_Mobster_file + " did not have a matching line in the input Mobster file " + args.in_Mobster + " and so this restart file must not be the correct restart file for this input Mobster file. Processing will not be carried out.\n")

	# Process each mobile-element-insertion-prediction record in the input Mobster file.
	# First process already-processed lines (if this is a restart) to be simply copied across to output, and header records.
	# Then process unprocessed lines, and process them in a linux process so that multiple processes (lines) can be carried out simultaneously on multi-core CPU.

	currently_processing_header_and_restart_lines = True
	restart_upto = 0
	in_upto = -1
	while (currently_processing_header_and_restart_lines):

		# Process header and already-processed restart lines

		in_upto = in_upto + 1
		if (in_upto >= len(input_Mobster_inlines)):
			currently_processing_header_and_restart_lines = False

		inline = input_Mobster_inlines[in_upto]
		wrote_out_this_record = False
		copy_this_line_across_from_restart_file = False

		if args.restart_Mobster_file is not None:
			if (in_upto < 4): # copy across header from restart file
				copy_this_line_across_from_restart_file = True
			else:
				# If there is a restart file and we haven't gotten to the end of it yet, then copy across restart lines instead of processing the matching input line,
				# because that input line was already processed by a previous run of this program and written out to the restart file, 
				# before that previous run crashed and didn't process all the input lines.
				# If there is a restart file and this input line is not in the restart file even though there are more lines in the restart file after where this line would be,
				# then this input line must have been dropped by the previous run of this program, due to the MEI being in an existing ME region and there being no evidence of an MEI in the BAM file.
				# In that case, we will process the input line again, unnecessarily, 
				# because it's easier/safer to do that than to try to match the input file to non-existant dropped lines in the restart file.
				if (restart_upto < len(restart_Mobster_inlines)):
					in_MEI_call = input_Mobster_inlines[in_upto]
					in_fields = in_MEI_call.split("\t")
					in_0 = in_fields[0]
					in_1 = in_fields[1]
					in_2 = in_fields[2]
					in_3 = in_fields[3]
					in_4 = in_fields[4]
					restart_MEI_call = restart_Mobster_inlines[restart_upto]
					restart_fields = restart_MEI_call.split("\t")
					restart_0 = restart_fields[0]
					restart_1 = restart_fields[1]
					restart_2 = restart_fields[2]
					restart_3 = restart_fields[3]
					restart_4 = restart_fields[4]	
					# If there exists a restart line that is the same as this input line, and it is a complete line,
					# then copy it across
					if ((in_0 == restart_0) and (in_1 == restart_1) and (in_2 == restart_2) and (in_3 == restart_3) and (in_4 == restart_4)):
						last_chr = restart_MEI_call[ len(restart_MEI_call)-1 : len(restart_MEI_call) ]
						last_chr_is_carriage_return = False
						if ((last_chr == "\n") or (last_chr == "\r")):
							last_chr_is_carriage_return = True
						if (last_chr_is_carriage_return):
							copy_this_line_across_from_restart_file = True
						else:
							restart_upto = restart_upto + 1 # we have processed this restart_file_line by deciding to ignore it, it isn't complete

		if (copy_this_line_across_from_restart_file):

			restart_inline = restart_Mobster_inlines[restart_upto]
			output_Mobster.write( restart_inline )
			restart_upto = restart_upto + 1

		else:

			inline = inline.replace( "\n", '' )
			inline = inline.replace( "\r", '' )
			if (inline != ''):
				if ((inline[0:1] == "#") or (inline[0:4] == "Chr\t")):

					# Write out header lines as-is
					# If it is the headings line, append field headings to it for the fields that this program may append to Mobster records

					outline = inline + "\n"
					if (inline[0:4] == "Chr\t"):
						outline = inline + "\t" + 'MEI_call_quality_from_Mobster_counts'
						outline = outline + "\t" + 'refined_border5' + "\t" + 'refined_border3'
						outline = outline + "\t" + 'refined_border5_MEI_sequence_fragment_start_pos' + "\t" + 'refined_border3_MEI_sequence_fragment_start_pos'
						outline = outline + "\t" + 'refined_border5_avg_read_count_for_MEI_read_fragment' + "\t" + 'refined_border3_avg_read_count_for_MEI_read_fragment'
						outline = outline + "\t" + 'refined_border5_length_of_MEI_read_fragment' + "\t" + 'refined_border3_length_of_MEI_read_fragment'
						outline = outline + "\t" + 'refined_border5_percent_that_non_MEI_part_of_MEI_fragments_matches_reference' + "\t" + 'refined_border3_percent_that_non_MEI_part_of_MEI_fragments_matches_reference'
						outline = outline + "\t" + 'refined_ref_seq_spanning_sequence_seen' + "\t" + 'refined_ref_seq_spanning_sequence_seen_start_pos'
						outline = outline + "\t" + 'refined_ref_seq_spanning_sequence_seen_end_pos' + "\t" + 'refined_ref_seq_spanning_sequence_seen_num_reads'
						outline = outline + "\t" + 'refined_border5_MEI_sequence_fragment' + "\t" + 'refined_border3_MEI_sequence_fragment' + "\n"
					output_Mobster.write( outline )
					wrote_out_this_record = True

					if args.discard_Mobster_file is not None:
						outline2 = inline + "\n"
						if (inline[0:4] == "Chr\t"):
							outline2 = 'Reason_for_discard' + "\t" + inline + "\t" + 'MEI_call_quality_from_Mobster_counts' + "\n"
						discard_Mobster.write( outline2 )

				else:

					# We are now up to processing mobile-element-insertion-prediction that need processing
					currently_processing_header_and_restart_lines = False

	currently_processing_new_prediction_lines = True
	if (in_upto >= len(input_Mobster_inlines)):
		currently_processing_new_prediction_lines = False
	while (currently_processing_new_prediction_lines):

		# We are now up to new mobile-element-insertion-predictions that need to be processed.

		num_concurrent_processes = 1
		if args.num_cores is not None:
			num_concurrent_processes = int(args.num_cores)
		quality_filter = ''
		if args.quality_filter is not None:
			quality_filter = args.quality_filter

		# Create an input batch of them to be processed by multiple concurrent processes.

		processing_batch_input = []
		for i in range( 0, num_concurrent_processes ):
			if ( (in_upto+i) < len(input_Mobster_inlines)):
				processing_batch_input.append( input_Mobster_inlines[(in_upto+i)] )

		# Send them off to be processed by multiple concurrent processes.

		processing_pool = Pool(num_concurrent_processes)
		processing_pool_result = processing_pool.map( process_one_input_line_from_processing_batch, processing_batch_input )
		processing_pool.close()

		# After multiple process processing of this input batch, write the outputs to the output files.

		for i in range( 0, len(processing_pool_result) ):
			this_result = processing_pool_result[i]
			if (this_result[0] != ''):
				output_Mobster.write( this_result[0] )
			if args.discard_Mobster_file is not None:
				if (this_result[1] != ''):
					discard_Mobster.write( this_result[1] )

		output_Mobster.flush()
		if args.discard_Mobster_file is not None:
			discard_Mobster.flush()

		in_upto = in_upto + num_concurrent_processes
		if (in_upto >= len(input_Mobster_inlines)):
			currently_processing_new_prediction_lines = False

	output_Mobster.close()
	if args.discard_Mobster_file is not None:
		discard_Mobster.close()


if __name__=='__main__':
    main()


