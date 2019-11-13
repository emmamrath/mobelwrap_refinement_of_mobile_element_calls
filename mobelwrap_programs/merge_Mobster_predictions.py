#!/usr/bin/python
# python merge_Mobster_predictions.py output_Mobster_predictions input_Mobster_1 input_Mobster_2 ...

# If Mobster predictions overlap, consider them to be the one prediction, remove duplicates.
# Choose the prediction that has the highest sample_counts.

# When is this program needed?
# When a Mobster job is too large to run (runs out of memory or wall-time)
# and so the input bam file to Mobster is split into multiple bam files and run separately.
# When split into multiple bam files by chromosome N, 
# must include the discordant reads on other chromosomes whose read pair is on chromosome N.

# Example of multiple predictions that must be merged into one prediction by choosing one of them only.
# Chr	Insert Point	border5		border3		sample_counts	cluster5	cluster3
# 5	144330992	144330978	144330994	test_sample=15	219		357
# 5	144330994	144330994	144331002	test_sample=43	470		357
# 5	144330994	144330994	144331013	test_sample=5	150		136

# Untested command to extract bam file reads by chromosome, including chromosome's discordant reads:
# awk '($1=="@") || (($3=="chr1") || ($7=="chr1") || ($3=="chr2") || ($7=="chr2"))'

# Example commands for a large sample to extract BAM file reads by chromosome,
# so that Mobster can be run on smaller BAM files
# whose results will be merged by this python program.
#
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1 ~ /^@/) || ((($3 == "1") || ($3 == "2") || ($3 == "3") || ($3 == "4")) || (($7 == "1") || ($7 == "2") || ($7 == "3") || ($7 == "4"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_chr1to4.bam
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1 ~ /^@/) || ((($3 == "5") || ($3 == "6") || ($3 == "7") || ($3 == "8")) || (($7 == "5") || ($7 == "6") || ($7 == "7") || ($7 == "8"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_chr5to8.bam
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1 ~ /^@/) || ((($3 == "9") || ($3 == "10") || ($3 == "11") || ($3 == "12")) || (($7 == "9") || ($7 == "10") || ($7 == "11") || ($7 == "12"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_chr9to12.bam
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1 ~ /^@/) || ((($3 == "13") || ($3 == "14") || ($3 == "15") || ($3 == "16")) || (($7 == "13") || ($7 == "14") || ($7 == "15") || ($7 == "16"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_chr13to16.bam
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1 ~ /^@/) || ((($3 == "17") || ($3 == "18") || ($3 == "19") || ($3 == "20")) || (($7 == "17") || ($7 == "18") || ($7 == "19") || ($7 == "20"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_chr17to20.bam
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1 ~ /^@/) || ((($3 == "21") || ($3 == "22") || ($3 == "X") || ($3 == "Y") || ($3 == "MT")) || (($7 == "21") || ($7 == "22") || ($7 == "X") || ($7 == "Y") || ($7 == "MT"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_chr21to22XYMT.bam
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam | awk -F'\t' '(($1~/^@/)||((($3!="1")&&($3!="2")&&($3!="3")&&($3!="4")&&($3!="5")&&($3!="6")&&($3!="7")&&($3!="8")&&($3!="9")&&($3!="10")&&($3!="11")&&($3!="12")&&($3!="13")&&($3!="14")&&($3!="15")&&($3!="16")&&($3!="17")&&($3!="18")&&($3!="19")&&($3!="20")&&($3!="21")&&($3!="22")&&($3!="X")&&($3!="Y")&&($3!="5")&&($3!="MT"))||(($7!="1")&&($7!="2")&&($7!="3")&&($7!="4")&&($7!="5")&&($7!="6")&&($7!="7")&&($7!="8")&&($7!="9")&&($7!="10")&&($7!="11")&&($7!="12")&&($7!="13")&&($7!="14")&&($7!="15")&&($7!="16")&&($7!="17")&&($7!="18")&&($7!="19")&&($7!="20")&&($7!="21")&&($7!="22")&&($7!="X")&&($7!="Y")&&($7!="MT"))))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH_noChrom.bam

# Here is an example command to remove hard clips from a BAM file.
# samtools view -h MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked.bam | awk -F'\t' '(($6 !~ /H/) || ($1 ~ /^@/))' | samtools view -Shu - > MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked_carefulRmvH.bam

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

# Mobster refers to:
# Genome Biol. 2014;15(10):488.
# Mobster: accurate detection of mobile element insertions in next generation sequencing data.
# Thung DT, de Ligt J, Vissers LE, Steehouwer M, Kroon M, de Vries P, Slagboom EP, Ye K, Veltman JA, Hehir-Kwa JY.

import sys
import os

######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def any_overlap( prediction_1, prediction_2 ):

	return_value = False

	prediction_1_strip = prediction_1.strip()
	infields_1 = prediction_1_strip.split("\t")

	prediction_2_strip = prediction_2.strip()
	infields_2 = prediction_2_strip.split("\t")

	if ((prediction_1_strip != '') and (prediction_2_strip != '')):

		prediction1_chrom = infields_1[0]
		prediction1_pos1 = infields_1[2]
		prediction1_pos2 = infields_1[3]
		prediction1_pos3 = infields_1[4]
		prediction2_chrom = infields_2[0]
		prediction2_pos1 = infields_2[2]
		prediction2_pos2 = infields_2[3]
		prediction2_pos3 = infields_2[4]

		if (prediction1_chrom == prediction2_chrom):
			start1 = min( prediction1_pos1, prediction1_pos2, prediction1_pos3 )
			end1 = max( prediction1_pos1, prediction1_pos2, prediction1_pos3 )
			start2 = min( prediction2_pos1, prediction2_pos2, prediction2_pos3 )
			end2 = max( prediction2_pos1, prediction2_pos2, prediction2_pos3 )
			if ((start1 <= start2) and (end1 >= start2)):
				return_value = True
			elif ((start1 <= start2) and (end1 >= end2)):
				return_value = True
			elif ((start2 <= start1) and (end2 >= start1)):
				return_value = True
			elif ((start2 <= start1) and (end2 >= end1)):
				return_value = True
	return return_value

######################################################
def best_sample( prediction_1, prediction_2 ):

	return_value = ''

	prediction_1_strip = prediction_1.strip()
	prediction_2_strip = prediction_2.strip()

	if ((prediction_1_strip != '') and (prediction_2_strip != '')):

		infields_1 = prediction_1_strip.split("\t")
		prediction1_sample_counts = infields_1[7]
		infields_1 = prediction1_sample_counts.split("=")
		prediction1_sample_counts = int(infields_1[1])

		infields_2 = prediction_2_strip.split("\t")
		prediction2_sample_counts = infields_2[7]
		infields_2 = prediction2_sample_counts.split("=")
		prediction1_sample_counts = int(infields_2[1])

		if ( prediction1_sample_counts >= prediction1_sample_counts ):
			return_value = prediction_1
		else:
			return_value = prediction_2

	else:
		if (prediction_2_strip != ''):
			return_value = prediction_2
		else:
			return_value = prediction_1
	return return_value

######################################################
def main():

	### read input files

	output_Mobster_file_name = sys.argv[1]
	input_Mobster_file_name = []
	if (len(sys.argv) > 2):
		for file_idx in range( 2, len(sys.argv) ):
			input_Mobster_file_name.append( sys.argv[file_idx] )

	all_input_Mobster_predictions_with_sort_key = []
	all_comment_lines = [] 		# Comment lines of all input files will be written to output file
	heading_line = '' 		# Heading line of only the first file will be written to output file

	for file_idx in range( 0, len(input_Mobster_file_name) ):

		input_Mobster = open(input_Mobster_file_name[file_idx], 'r')
		Mobster_lines = input_Mobster.readlines()
		input_Mobster.close()
		for line_idx in range( 0, len(Mobster_lines) ):
			inline = Mobster_lines[line_idx]
			is_comment_line = False
			is_heading_line = False
			strip_inline = inline.strip()
			if (strip_inline != ''):
				if (len(strip_inline) > 3):
					if (strip_inline[0:1] == '#'):
						all_comment_lines.append( inline )
					elif (strip_inline[0:4] == "Chr\t"):
						if (heading_line == ''):
							heading_line = inline
					else:

						infields = strip_inline.split("\t")
						chrom = infields[0]
						pos1 = infields[2]
						pos2 = infields[3]
						pos3 = infields[4]

						key_chrom = "%20s" % chrom
						key_pos1 = "%20d" % int(pos1)
						key_pos2 = "%20d" % int(pos2)
						key_pos3 = "%20d" % int(pos3)
						key_chrom = key_chrom[0:20]
						key_pos1 = key_pos1[0:20]
						key_pos2 = key_pos2[0:20]
						key_pos3 = key_pos3[0:20]
						sort_key = key_chrom + key_pos1 + key_pos2 + key_pos3
						inline_with_key = sort_key + inline
						all_input_Mobster_predictions_with_sort_key.append( inline_with_key )

	# sort and choose one unique prediction for each group of one or more overlapping predictions

	all_input_Mobster_predictions_with_sort_key.sort()

	all_input_Mobster_predictions = []
	for prediction_with_sort_key in all_input_Mobster_predictions_with_sort_key:
		prediction = prediction_with_sort_key[80:]
		all_input_Mobster_predictions.append( prediction )

	best_unique_Mobster_predictions = []

	previous_best_inline = ''
	for line_idx in range( 0, len(all_input_Mobster_predictions) ):
		inline = all_input_Mobster_predictions[line_idx]
		any_overlap_result = any_overlap( inline, previous_best_inline )
		if (any_overlap_result): # in a group of overlapping predictions, choose the best one
			previous_best_inline = best_sample( inline, previous_best_inline )
		else: # not yet or no longer in a group of overlapping predictions
			if (previous_best_inline != ''):
				best_unique_Mobster_predictions.append( previous_best_inline )
			previous_best_inline = inline
	if (previous_best_inline != ''):
		best_unique_Mobster_predictions.append( previous_best_inline )

	# output the list of unique predictions, having chosen the best one of any duplicates

	output_Mobster = open(output_Mobster_file_name, 'w')

	for line_idx in range( 0, len(all_comment_lines) ):
		outline = all_comment_lines[line_idx]
		output_Mobster.write( outline )
		
	if (heading_line != ''):
		output_Mobster.write( heading_line )

	for line_idx in range( 0, len(best_unique_Mobster_predictions) ):
		outline = best_unique_Mobster_predictions[line_idx]
		output_Mobster.write( outline )

	output_Mobster.close()


if __name__=='__main__':
    main()


