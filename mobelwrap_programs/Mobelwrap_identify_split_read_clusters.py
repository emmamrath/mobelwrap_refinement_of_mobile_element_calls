#!/usr/bin/python
# python Mobelwrap_identify_split_read_clusters.py what_is_the_extract where_is_bam [read_length] [indexed_reference_genome_fasta_file]
# python Mobelwrap_identify_split_read_clusters.py 1:564300-564400 temp_MGRBp1_sample1000_1_564300_564400.bam > identify_sample1000_1_564300_564400.txt
# python Mobelwrap_identify_split_read_clusters.py 1:564300-564400 temp_MGRBp1_sample1000_1_564300_564400.bam 150 > identify_sample1000_1_564300_564400.txt
# python Mobelwrap_identify_split_read_clusters.py 1:564300-564400 temp_MGRBp1_sample1000_1_564300_564400.bam - /home/biomath/projects_associated_with_mobile_elements_and_telomeres_and_mitochondrial_DNA/fasta_of_hg19/human_g1k_v37.fasta > identify_sample1000_1_564300_564400.txt
# python Mobelwrap_identify_split_read_clusters.py 1:564336-564704 /nvme/mobile_elements_2016dec/debug_MGRB_picard_header/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.bam 150

# In a small region of the genome, in a bam file of reads mapped to that region,
# we may see that all reads map nicely to the reference genome, regardless of the read start position.
# Or else, we may see they some of the reads map nicely for the entire read
# while others map nicely for only a part of the read and the rest of the read does not map here.
# These reads will have been mapped with a soft clip in their bam file cigar code.
# This situation occurs for mobile element insertions and for chromosomal rearrangements.
# For a mobile element insertion, in the bam file we see a continuous ref. seq. of read
# and also two other sequences for the 2 ends of the broken ref. seq. having mobile element sequence attached to it.

# This program takes a region as input, and identifies the different continuous sequences seen in the bam file
# when stringing overlapping reads together.
# This program writes its output to STDOUT.
# It writes the chromosome and starting_position on separate lines.
# Then it writes out the list of unique concatenated bam sequences it has seen, one per line.
# The starting_position that this program writes out is usually different to the input starting position
# because the program receives from samtools and inspects bam reads that overlap the input starting position yet may extend beyond the input positions.

# The indexed_reference_genome_fasta_file needs to have been indexed with samtools faidx.

# For an input region of chrom1:pos1-pos2, this program looks at the region chrom1:(pos1-read_length):(pos1+read_length).

# If read_length is not given on input, then this program first assumes that it is 150.
# If read_length is not given on input and after reading bam records this program finds that the bam reads are longer, 
# this program will set the read_length to the larger value, and re-read bam records from a larger region.

# When using this program to identify unique concatenated reads in an area where there is a mobile element insertion,
# the smaller input region to this program, the better.
# Longer input regions to this program increase the occurrance of the program reading bam reads that have a 1 bp deletion compared to the other reads,
# and this manifests itself as a unique different concatenated sequence.
# Really it is not a new sequence. 

# This program tries to handle small shifts caused by 1+ bp INDELs.
# It compares parts of the read to the reference genome.
# If they match the reference genome at a position that is only 1 to 5 off the reference genome,
# then the read is discarded as a stutter.

# This program recognises some stutters of mobile-element-insertions not matching the reference genome
# only if caused by poly-A or poly-T sequences.
# Otherwise, this program does not recognise stutters of mobile-element-insertions not matching the reference genome
# and they will be identified as multiple unique sequences, each having only 1 read record (each one is a weak MEI candidate).
# What should happen is that they be recognised as the same mobile-element-insertion sequence having multiple reads supporting it (makes it a strong MEI candidate).

# For example, here is an mobile element insertion at 1:13410-813915:
#
#    1   11830078    143M7S      2   234632297     CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTT
#    1   11830080    141M9S      10  101597212     CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTT
#    1   11830080    141M9S      19  30389229      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTT
#    1   11830082    139M10S     =   11830082      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTT
#    1   11830082    *           =   11830082      TCCGAGTGCGCGGGGATAGTTCTTAAGGATACGGCGATCGCCGAGATCTCCCCCAGGCCGTACACAATAAACAAAAAAGAAGCCCTCCCGAT
#    1   11830083    138M12S     =   11829918      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTT
#    1   11830085    149M        =   11830333      CCCTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAA
#    1   11830089    150M        =   11830370      CCCTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGT
#    1   11830089    150M        =   11830370      CCCTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGT
#    1   11830096    125M25S     =   11829651      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTAAATTTATTT
#    1   11830100    121M29S     17  20889245      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTTAATTTATTTTTTT
#    1   11830106    115M35S     10  101597212     CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTAAATTTATTTTTTTTTTGAT
#    1   11830109    112M38S     19  30389452      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTAAATTTATTTTTTTATTGATAAT
#    1   11830109    112M38S     12  104716787     CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTAAATTTATTTTTTTATTGATAAT
#    1   11830111    110M39S     =   11829960      CCCTGGAAGTTTTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTTTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTAAATTTATTTTTTTATTGATAATTC
#    1   11830112    109M39S     =   11829960      CCCTGGAAGTATTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTAAATTTATTTTTTTATTGATAATTC
#    1   11830117    104M46S     11  125734030     CCCTGGAAGTCTTTTTATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATCTTAAATTTTTCTTTTTTCTTTTTTTTTTATACTTATTGTTGTCTGGGTAATGTTTTCGT
#    1   11830126    95M55S      2   234632270     CCCTGGAAGTATTTTTATGAATATTTGCACCGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTAAATCTATTTTTTTATGGATAATTCTTGGGCGTTTCTCAC
#    1   11830133    150M        =   11829727      CCTTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACA
#    1   11830144    147M        =   11830144          GGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830144    *           =   11830144          ACCCGGGGTCGCCCGCACCCCACGTAACCTGCCGACTGGGGGCGATGACCCCCCCACCTCGTCGGAGACCACCAACACAACAAAAAAAGAAAACAACCAGCCAGGGGGCCGGGGGCGGTCATCTCGCGCCCAAGGTCCTCGAGAGC
#    1   11830151    149M        =   11830385                 TTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830151    149M        =   11830385                 TTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830151    149M        =   11830385                 TTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830152    69M80H      =   11829592                TTTTATGAATGTTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTT
#    1   11830154    149M        =   11830369                    TTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830155    63M84H      17  64265621                     TATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTT
#    1   11830155    66M83H      17  64265621                     TATGAATATTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTT
#    1   11830163    58M91H      19  30389391                             TTTGCACAGCATTTTTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTT
#    1   11830170    37S23M90S   19  22987044      AAAGGGGGCAATCTCCGGAATTTTTAAAGCAGCATTTTTTAAAAAAAAGGGAATTTGGCACCTCCTTTTTTCACACATCGGTTTTCATAAGAACTTTCTTTGTGTTTTTCATTCATCCCTTTATAGAAATTTACAAAAGGATG
#    1   11830173    150M        =   11830674                                       ATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830173    150M        =   11830402                                      ATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGCGCA
#    1   11830175    149M        =   11830344                                         TTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#                                                  0      0150      0160      0170      0180      0190      0200      0210      0220      0230      0240      0250      0260      0270      0280      029
#                                                  4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----
#                                                  CCCTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830175    149M        =   11830344                                         TTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830177    44M104H     15  35949380                                           TTTTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTT
#    1   11830177    149M        =   11830011                                           TTTAAAAAAAAGGGAAAATGGTGTCGCAATTTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTCCCCACCCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830177    149M        =   11830011                                           TTTAAAAAAAAGGGAAAATTGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830179    42M108H     12  104716795                                            TTAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTT
#    1   11830181    150M        =   11830537                                               AAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830188    150M        =   11830677                                                     GGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAAGATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830189    150M        =   11829900                                                       GGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCCCAAAGTTGCATAGAAACCACATGGTGCA
#    1   11830196    150M        =   11830452                                                              GGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830198    38M1D25M86S =   11830533                                                                TGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAAGTAACCAAAATAATTTTCCCAACCAAATAAGTTCCAAAGAAACAAAAGGGGGCT
#    1   11830209    22S128M     =   11830570                                                     CGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    6S144M      =   11830543                                                                     GCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    56S94M      =   11830348                   CTGCTCCTTGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGGA
#    1   11830209    79H71M      =   11830344                                                                           TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACC
#    1   11830209    100H50M     8   9515677                                                                           TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAA
#    1   11830209    88H62M      =   11830649                                                                           TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGC
#    1   11830209    82H67M      =   11830442                                                                          TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGA
#    1   11830209    82H68M      9   91240906                                                                           TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAA
#    1   11830209    100H50M     9   115276084                                                                          TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAA
#    1   11830209    113H37M     17  15540391                                                                           TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAA
#    1   11830209    105H45M     hs37d56866176                                                                          TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTT
#    1   11830209    62S87M      =   11830314             AGGCGGCTGCTCCTTGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    62S87M      =   11830314             AGGCGGCTGCTCCTTGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    48S101M     16  56774327                           TGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    81H68M      =   11830386                                                                          TAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACCTAGTTGCATACAA
#    1   11830209    48S100M     15  35947991                           TGCCACCGGGCCCCGCGGGGCCCGCCCGCCCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    48S100M     15  35947993                           TGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    44S105M     =   11830525                             CTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    44S105M     =   11830525                               CTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830209    48S101M     =   11830365                          TGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830214    150M        =   11830484                                                                                TTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830217    149M        =   11830509                                                                                TCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830219    150M        =   11830563                                                                                   TTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830220    148M        18  46728812                                                                                     TAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCCCATAGTTGCATAGAAACCACAGGGTGCA
#    1   11830220    148M        11  119379015                                                                                    TAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830221    149M        =   11830531                                                                                    AAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830226    150M        8   65597807                                                                                        AGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830229    150M        =   11830571                                                                                              ATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830245    150M        =   11830542                                                                                                         CCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#    1   11830248    150M        =   11830517                                                                                                             AATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
 
# By visually joining overlapping reads, we see 3 different sequences, 
# which are a mobile element (ME) insertion 3'end (polyT tail), the reference sequence, and the ME insertion 5'end:
#
# CCCTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTTTTTCTTTTTTTTTTAAATTTATTTTTTTATTGATAAT
# CCCTGGAAGTATTTTTATGAATATTTGCATAGCATTTTTTAAAAAAAAGGGAAAATGGTGTGGCAATGTTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA
#        AGGCGGCTGCTCCTTGCCCTCGGGCCCCGCGGGGCCCGTCCGCTCCTCCAGCCGCTGCCTCCTAAATTTTTCTTAAAGGAGCATAAAAATGTATCCAAAATATTTTTCCCAACCACATAGTTGCATAGAAACCACATGGTGCA

# The goal of this program is to produce those 3 different sequences seen in the input.

# This program uses samtools and command line wc -l to find out how many sequences in a region, eg:
# samtools view MGRBp1_sample630_sg1_humansingle1_874.picard_sorted.dupmarked.bam "2:31951000-31953000" | wc -l
# If there are too many reads in a region, then the MEI is not inspected, 
# and the IMPRECISE Mobster MEI prediction for which this program is being called will remain IMPRECISE, won't because PRECISE from inspection by this program.

# Some of the functions in this program has been broken into smaller functions so that they can be profiled with:
# python -m cProfile -s cumulative identify_split_read_clusters_SPLIT_INTO_MORE_FUNCTIONS.py 1:566220-566280 \
#	temp_MGRBp1_sample1000_1_566262_813815.bam - \
#	/home/biomath/projects_associated_with_mobile_elements_and_telomeres_and_mitochondrial_DNA/reference_genomes/hs37d5x/hs37d5x.fa \
#	> sample1000_1_566220_566280_identifyClustersFUNCTIONS_output.txt > sample1000_1_566220_566280_identifyClustersFUNCTIONS_cProfile_output.txt



__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'


import sys
import os
import datetime
import math
import subprocess
import commands
import argparse
import re


ref_seq_chrom = ''
ref_seq_start_pos = 0
ref_seq_length = 0
ref_seq_array = []
bam_read_chrom = []
bam_read_pos = []
bam_read_cigar = []
bam_read_sequence = []
bam_read_num_reads = 0
# bam_stutter_read_chrom = []
# bam_stutter_read_pos = []
# bam_stutter_read_cigar = []
# bam_stutter_read_sequence = []
bam_read_longest_read_length = 0
seen_sequence = []
seen_sequence_count = []
final_seen_sequence = []
final_seen_sequence_count = []

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
def collapse_seq_count( in_seq ):

	# The result will be passed back to another program (convert_imprecise_Mobster_predictions_into_precise.py)
	# Each base pair must occupy 1 and only 1 character.
	# So turn counts > 9 into a '9'

	out_seq = ''
	for this_bp in in_seq:
		if (this_bp > 9):
			out_seq = out_seq + '9'
		else:
			this_bp = str(this_bp)
			this_bp = this_bp[0:1]
			out_seq = out_seq + this_bp
	return out_seq

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
def is_this_read_a_stutter_of_ref_seq( read_cigar ):

	# is it a pattern such as 14M1D136M
	is_a_stutter = False
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
			is_a_stutter = True

	# is it a pattern such as 100M1I48M
	if (is_a_stutter == False):
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
				is_a_stutter = True

	return is_a_stutter

######################################################
def get_bam_reads( where_is_bam, chrom, num1, num2 ):

	global bam_read_chrom
	global bam_read_pos
	global bam_read_cigar
	global bam_read_sequence
	global bam_read_num_reads
	global bam_read_longest_read_length

	# How many reads are in this MEI region of the genome?
	# If there are too many, then it is probably a repetitive-sequence region
	# and trying to figure out the MEI positions will take a long time and be inaccurate.
	# We seem to be seeing with the Illumina HighTenX that such regions have 100x more reads than other regions.
	# So let's not process MEIs in regions having lots of reads.

        max_reads_per_bp = 2
        count_reads_command = 'samtools view -c ' + where_is_bam + ' ' + str(chrom) + ':' + str(num1) + '-' + str(num2)
	count_command_status, count_command_output = commands.getstatusoutput( count_reads_command )
	if (count_command_status != 0):
		raise ValueError("\n\nIn Mobelwrap_identify_split_read_clusters.py, command failed with status: " + str(count_command_status) + " for command: " + count_reads_command + "\nThus will not continue processing.\n")
	else:
		count_command_output_lines = parse_out_any_warnings( count_command_output, count_reads_command )
		if (is_integer(count_command_output_lines[0]) == False):
			raise ValueError("\n\nIn Mobelwrap_identify_split_read_clusters.py, command failed with output: " + str(count_command_output) + " for command: " + count_reads_command + "\nThus will not continue processing.\n")
		else:
			# This line contains an integer and is our result
			count_reads_output = int(count_command_output_lines[0])
			actual_reads_per_bp = float(count_reads_output) / float(num2 - num1 + 1)

        if (actual_reads_per_bp <= max_reads_per_bp):

		# get all the reads in this MEI region of the genome

		outlines = []
		samtools_command = 'samtools view ' + where_is_bam + ' ' + str(chrom) + ':' + str(num1) + '-' + str(num2)
		samtools_status, samtools_output = commands.getstatusoutput( samtools_command )
		if (samtools_status != 0):
			raise ValueError("\n\nIn Mobelwrap_identify_split_read_clusters.py, command failed with status: " + str(samtools_status) + " for command: " + samtools_command + "\nThus will not continue processing.\n")
		else:
			outlines = parse_out_any_warnings( samtools_output, samtools_command )

		bam_read_chrom = [''] * len(outlines)
		bam_read_pos = [-1] * len(outlines)
		bam_read_cigar = [''] * len(outlines)
		bam_read_sequence = [''] * len(outlines)
		bam_upto = 0
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

				# adjust the position of this bam read for its soft clips

				first_soft_clip = 0
				cigar_split = read_cigar.split('S')
				if (is_integer(cigar_split[0])):
					first_soft_clip = int(cigar_split[0])
				read_pos = read_pos - first_soft_clip

				keep_this_sequence = True
				if (read_cigar == '*'):	# sequences that don't align to reference sequence will have cigar equals *
					keep_this_sequence = False

				if (read_cigar.find('H') != -1): # Don't look at hard clips. There will be a soft clip version of this read that we will look at.
					keep_this_sequence = False

				if (keep_this_sequence):

					# discard sequences that are a stutter of the reference sequence

					is_it_a_stutter = is_this_read_a_stutter_of_ref_seq( read_cigar )
					if (is_it_a_stutter == False):
						bam_read_chrom[ bam_upto ] = read_chrom
						bam_read_pos[ bam_upto ] = read_pos
						bam_read_cigar[ bam_upto ] = read_cigar
						bam_read_sequence[ bam_upto ] = read_sequence
						bam_upto = bam_upto + 1
		bam_read_num_reads = bam_upto

	return

######################################################
def are_these_2_sequences_the_same_PREPARE_COMPARISON_STRINGS( seq1, seq2 ):

	match_percent = 0
	match_percent_5prime = 0
	match_percent_3prime = 0

	record_matches_mismatches = [' '] * len(seq1)
	for i in range( 0, len(seq1) ):
		if ((seq1[i] != ' ') and (seq2[i] != ' ')):
			if (seq1[i] == seq2[i]):
				record_matches_mismatches[i] = '1'
			else:
				record_matches_mismatches[i] = '0'

	collapse_seq1 = collapse_seq(seq1)
	collapse_seq2 = collapse_seq(seq2)
	collapse_matches_mismatches = collapse_seq(record_matches_mismatches)
	collapse_matches_mismatches = collapse_matches_mismatches.replace(' ','')
	count_matches = len(collapse_matches_mismatches.replace('0',''))
	count_bp = len(collapse_matches_mismatches)
	if (count_bp == 0):
		match_percent = 0
		match_percent_5prime = 0
		match_percent_3prime = 0
	else:
		match_percent = count_matches * 100 / count_bp
		length_of_each_end_to_compare = 20
		match_percent_5prime_string = collapse_matches_mismatches[0:length_of_each_end_to_compare]
		end1 = len(collapse_matches_mismatches) - length_of_each_end_to_compare
		end2 = len(collapse_matches_mismatches)
		if (end1 < 0):
			end1 = 0
		if (end2 < 0):
			end2 = 0
		match_percent_3prime_string = collapse_matches_mismatches[end1:end2]
		count_matches_5prime = len(match_percent_5prime_string.replace('0',''))
		count_matches_3prime = len(match_percent_3prime_string.replace('0',''))
		if (len(match_percent_5prime_string) == 0):
			match_percent_5prime = 0
		else:
			match_percent_5prime = count_matches_5prime * 100 / len(match_percent_5prime_string)
		if (len(match_percent_3prime_string) == 0):
			match_percent_3prime = 0
		else:
			match_percent_3prime = count_matches_3prime * 100 / len(match_percent_3prime_string)

	return record_matches_mismatches, collapse_seq1, collapse_seq2, match_percent, match_percent_5prime, match_percent_3prime

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_STUTTERS( sameness_result, AAAAAA_or_TTTTTT, collapse_seq1, collapse_seq2, minimum_match_percent ):

	# These 2 sequences don't match exactly.
	# Let's see if they are stutter sequences of each other that don't match exactly due to poly-A or poly-T sequences of a mobile-element-insertion

	match_percent = 0
	AAAAAA_or_TTTTTT_plus = AAAAAA_or_TTTTTT + '+'

	seq1_collapse_to_AAAAAA = re.sub( AAAAAA_or_TTTTTT_plus, AAAAAA_or_TTTTTT, collapse_seq1 )
	seq2_collapse_to_AAAAAA = re.sub( AAAAAA_or_TTTTTT_plus, AAAAAA_or_TTTTTT, collapse_seq2 )
	seq1_collapse_to_AAAAAA = seq1_collapse_to_AAAAAA.strip()
	seq2_collapse_to_AAAAAA = seq2_collapse_to_AAAAAA.strip()
	seq1_find_AAAAAA = seq1_collapse_to_AAAAAA.find( AAAAAA_or_TTTTTT )
	seq2_find_AAAAAA = seq2_collapse_to_AAAAAA.find( AAAAAA_or_TTTTTT )
	if ((seq1_find_AAAAAA >= 0) and (seq2_find_AAAAAA >= 0)):
		difference_in_length = abs( seq1_find_AAAAAA - seq2_find_AAAAAA )
		if (seq1_find_AAAAAA > seq2_find_AAAAAA):
			seq2_collapse_to_AAAAAA = (' ' * difference_in_length) + seq2_collapse_to_AAAAAA
		elif (seq2_find_AAAAAA > seq1_find_AAAAAA):
			seq1_collapse_to_AAAAAA = (' ' * difference_in_length) + seq1_collapse_to_AAAAAA
		longest_length = len(seq1_collapse_to_AAAAAA)
		if (len(seq2_collapse_to_AAAAAA) < longest_length):
			difference_in_length = longest_length - len(seq2_collapse_to_AAAAAA)
			seq2_collapse_to_AAAAAA = seq2_collapse_to_AAAAAA + (' ' * difference_in_length)
		elif (len(seq2_collapse_to_AAAAAA) > longest_length):
			longest_length = len(seq2_collapse_to_AAAAAA)
			difference_in_length = longest_length - len(seq1_collapse_to_AAAAAA)
			seq1_collapse_to_AAAAAA = seq1_collapse_to_AAAAAA + (' ' * difference_in_length)
		count_matches = 0
		count_bp = 0
		for i in range( 0, longest_length ):
			count_bp = count_bp + 1
			if ( seq1_collapse_to_AAAAAA[i:i+1] == seq2_collapse_to_AAAAAA[i:i+1] ):
				count_matches = count_matches + 1
		match_percent = count_matches * 100 / count_bp
		if (match_percent >= minimum_match_percent):
			sameness_result = 'SAME_BUT_WITH_STUTTER'

	return sameness_result, match_percent

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_5PRIME_MATCH( sameness_result, record_matches_mismatches, minimum_number_of_matches ):

	# These 2 sequences don't match exactly.
	# Let's see if they match on the 5' end or on the 3' end.
	# If so, then it is a split read with part of the sequence matching.
	# Try to find the right-most position at which 90% of positions match on the left.
	# Try to find the left-most position at which 90% of positions match on the right.

	looking_for_5prime_match = True
	count_matches = 0
	# count_mismatches = 0
	total_count = 0
	the_3prime_edge_of_5prime_match = -1
	for i in range( 0, len(record_matches_mismatches) ):
		if (looking_for_5prime_match):
			if (record_matches_mismatches[i] != ' '):
				total_count = total_count + 1
				if (record_matches_mismatches[i] == '1'):
					count_matches = count_matches + 1
				# else:
					# count_mismatches = count_mismatches + 1
				if (total_count > minimum_number_of_matches):
					if ( (float(count_matches)/float(total_count)) >= 0.9):
						# looking_for_5prime_match is still true
						the_3prime_edge_of_5prime_match = i
					else: # we have seen too many non-matches, so there is no 5prime matching region
						looking_for_5prime_match = False
	if (the_3prime_edge_of_5prime_match > len(record_matches_mismatches)-4):
		sameness_result = 'SAME'

	return sameness_result, the_3prime_edge_of_5prime_match

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_3PRIME_MATCH( sameness_result, record_matches_mismatches, minimum_number_of_matches, minimum_match_fraction ):

	# These 2 sequences don't match exactly.
	# Let's see if they match on the 5' end or on the 3' end.
	# If so, then it is a split read with part of the sequence matching.
	# Try to find the right-most position at which 90% of positions match on the left.
	# Try to find the left-most position at which 90% of positions match on the right.

	looking_for_3prime_match = True
	count_matches = 0
	# count_mismatches = 0
	total_count = 0
	the_5prime_edge_of_3prime_match = len(record_matches_mismatches)
	for i in range( len(record_matches_mismatches)-1, -1, -1 ):
		if (looking_for_3prime_match):
			if (record_matches_mismatches[i] != ' '):
				total_count = total_count + 1
				if (record_matches_mismatches[i] == '1'):
					count_matches = count_matches + 1
				# else:
					# count_mismatches = count_mismatches + 1
				if (total_count > minimum_number_of_matches):
					compare_length = float(count_matches)/float(total_count)
					compare_last10_count = 0
					for j in range( i, i+9 ):
						if (i < len(record_matches_mismatches)):
							if (record_matches_mismatches[j] == '1'):
								compare_last10_count = compare_last10_count + 1
					compare_last10 = float(compare_last10_count)/float(10)
					compare_last1 = float(0)
					if (record_matches_mismatches[i] == '1'):
						compare_last1 = float(1)
					if ((compare_length >= minimum_match_fraction) and (compare_last10 >= minimum_match_fraction)):
						if (compare_last1 == 1):
							# looking_for_3prime_match is still true
							the_5prime_edge_of_3prime_match = i
					else: # we have seen too many non-matches, so there is no 5prime matching region
						looking_for_3prime_match = False
	if (the_5prime_edge_of_3prime_match < 3):
		sameness_result = 'SAME'

	return sameness_result, the_5prime_edge_of_3prime_match

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE_tally_matches_at_each_pos( record_matches_mismatches ):

	pos_having_highest_matches = -1

	int_record_matches_mismatches = [0] * len(record_matches_mismatches)
	for i in range( 0, len(record_matches_mismatches) ):
		if (record_matches_mismatches[i] == '1'):
			int_record_matches_mismatches[i] = 1
	num_matches_at_each_pos = [0] * len(record_matches_mismatches)

	for i in range( 0, 11 ):
		pos_to_add = i + 10
		num_matches_at_each_pos[i] = num_matches_at_each_pos[i-1] + int_record_matches_mismatches[pos_to_add]
	for i in range( 11, (len(record_matches_mismatches) - 10) ):
		fill_pos = i
		prev_pos_to_use = i - 1
		pos_to_drop = i - 11
		pos_to_add = i + 10
		num_matches_at_each_pos[fill_pos] = num_matches_at_each_pos[prev_pos_to_use] - int_record_matches_mismatches[pos_to_drop] + int_record_matches_mismatches[pos_to_add]
	for i in range( (len(record_matches_mismatches) - 10), len(record_matches_mismatches) ):
		pos_to_drop = i - 11
		num_matches_at_each_pos[i] = num_matches_at_each_pos[i-1] - int_record_matches_mismatches[pos_to_drop]

	pos_having_highest_matches = -1
	highest_matches = -1
	for i in range( 0, len(num_matches_at_each_pos) ):
		if (num_matches_at_each_pos[i] > highest_matches):
			pos_having_highest_matches = i
			highest_matches = num_matches_at_each_pos[i]

	return pos_having_highest_matches

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE_for_3prime_end( record_matches_mismatches, pos_having_highest_matches, the_3prime_end ):

	keep_looking_for_3prime_end = True
	count_matches_on_3prime_side = 0
	for i in range( pos_having_highest_matches, -1, -1 ):
		if (keep_looking_for_3prime_end == True):
			if (record_matches_mismatches[i] == '1'):
				count_matches_on_3prime_side = count_matches_on_3prime_side + 1
				the_3prime_end = i
			else:
				keep_looking_for_3prime_end = False

	return the_3prime_end

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE_for_5prime_end( record_matches_mismatches, pos_having_highest_matches, the_5prime_end ):

	keep_looking_for_5prime_end = True
	count_matches_on_5prime_side = 0
	for i in range( pos_having_highest_matches, len(record_matches_mismatches) ):
		if (keep_looking_for_5prime_end == True):
			if (record_matches_mismatches[i] == '1'):
				count_matches_on_5prime_side = count_matches_on_5prime_side + 1
				the_5prime_end = i
			else:
				keep_looking_for_5prime_end = False

	return the_5prime_end

######################################################
def are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE( sameness_result, sameness_result_pos, sameness_result_pos2, record_matches_mismatches ):

	# Let's see if the two sequences match in the middle.
	# For each position, record how many matches when look at 10 bp on either side.
	# If this position is the centre of at least 21 bp having 100% match, then it is a match in the middle.

	the_5prime_end = -1
	the_3prime_end = -1

	pos_having_highest_matches = are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE_tally_matches_at_each_pos( record_matches_mismatches )

	if (pos_having_highest_matches > -1):
		the_3prime_end = len(record_matches_mismatches)
		the_5prime_end = -1

		the_3prime_end = are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE_for_3prime_end( record_matches_mismatches, pos_having_highest_matches, the_3prime_end )

		the_5prime_end = are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE_for_5prime_end( record_matches_mismatches, pos_having_highest_matches, the_5prime_end )

		if ( ((the_5prime_end - the_3prime_end + 1) >= 21) and (the_3prime_end > 0) ):
			sameness_result = 'MATCHING_MIDDLE'
			sameness_result_pos = the_5prime_end
			sameness_result_pos2 = the_3prime_end

	return sameness_result, sameness_result_pos, sameness_result_pos2, the_5prime_end, the_3prime_end

######################################################
def are_these_2_sequences_the_same( seq1, seq2 ):

	sameness_result = 'DIFFERENT'
	sameness_result_pos = ''
	sameness_result_pos2 = ''
	minimum_number_of_matches = 10
	minimum_match_percent = 90
	minimum_match_percent_at_ends = 95
	minimum_match_fraction = 0.9

	record_matches_mismatches, collapse_seq1, collapse_seq2, match_percent, match_percent_5prime, match_percent_3prime = are_these_2_sequences_the_same_PREPARE_COMPARISON_STRINGS( seq1, seq2 )

	if ((match_percent >= minimum_match_percent) and (match_percent_5prime >= minimum_match_percent_at_ends) and (match_percent_3prime >= minimum_match_percent_at_ends)):

		# These 2 sequences match exactly or almost exactly along the length of positions for which both sequences have a nucleotide present.

		sameness_result = 'SAME'

	if (sameness_result != 'SAME'):

		# These 2 sequences don't match exactly.
		# Let's see if they are stutter sequences of each other that don't match exactly due to poly-A or poly-T sequences of a mobile-element-insertion

		sameness_result, match_percent = are_these_2_sequences_the_same_LOOK_FOR_STUTTERS( sameness_result, 'AAAAAA', collapse_seq1, collapse_seq2, minimum_match_percent )

	if ((sameness_result != 'SAME') and (sameness_result != 'SAME_BUT_WITH_STUTTER')):

		# These 2 sequences don't match exactly.
		# Let's see if they are stutter sequences of each other that don't match exactly due to poly-A or poly-T sequences of a mobile-element-insertion

		sameness_result, match_percent = are_these_2_sequences_the_same_LOOK_FOR_STUTTERS( sameness_result, 'TTTTTT', collapse_seq1, collapse_seq2, minimum_match_percent )

	if ((sameness_result != 'SAME') and (sameness_result != 'SAME_BUT_WITH_STUTTER')):

		# These 2 sequences don't match exactly.
		# Let's see if they match on the 5' end or on the 3' end.
		# If so, then it is a split read with part of the sequence matching.
		# Try to find the right-most position at which 90% of positions match on the left.
		# Try to find the left-most position at which 90% of positions match on the right.

		sameness_result, the_3prime_edge_of_5prime_match = are_these_2_sequences_the_same_LOOK_FOR_5PRIME_MATCH( sameness_result, record_matches_mismatches, minimum_number_of_matches )

		sameness_result, the_5prime_edge_of_3prime_match = are_these_2_sequences_the_same_LOOK_FOR_3PRIME_MATCH( sameness_result, record_matches_mismatches, minimum_number_of_matches, minimum_match_fraction )

		if (sameness_result == 'DIFFERENT'):

			# Let's see if the two sequences match in the middle.
			# For each position, record how many matches when look at 10 bp on either side.
			# If this position is the centre of at least 21 bp having 100% match, then it is a match in the middle.

			sameness_result, sameness_result_pos, sameness_result_pos2, the_5prime_end, the_3prime_end = are_these_2_sequences_the_same_LOOK_FOR_MATCHING_MIDDLE( sameness_result, sameness_result_pos, sameness_result_pos2, record_matches_mismatches )

		if (sameness_result == 'DIFFERENT'):

			if ((the_3prime_edge_of_5prime_match > -1) and (the_5prime_edge_of_3prime_match < len(record_matches_mismatches))):
				if ( (the_3prime_edge_of_5prime_match < the_5prime_edge_of_3prime_match) and ((the_5prime_edge_of_3prime_match - the_3prime_edge_of_5prime_match + 1) >= 21) ):
					sameness_result = 'MATCHING_EDGES_NOT_MIDDLE'
					sameness_result_pos = the_5prime_edge_of_3prime_match
					sameness_result_pos2 = the_3prime_edge_of_5prime_match
				elif ( (the_3prime_edge_of_5prime_match > the_5prime_edge_of_3prime_match) and ((the_3prime_edge_of_5prime_match - the_5prime_edge_of_3prime_match + 1) >= 21) ):
					sameness_result = 'MATCHING_MIDDLE'
					sameness_result_pos = the_5prime_edge_of_3prime_match
					sameness_result_pos2 = the_3prime_edge_of_5prime_match

		if (sameness_result == 'DIFFERENT'):

			if (the_3prime_edge_of_5prime_match > -1):
				sameness_result = '5_PRIME_SAME'
				sameness_result_pos = the_3prime_edge_of_5prime_match
			elif (the_5prime_edge_of_3prime_match < len(record_matches_mismatches)):
				sameness_result = '3_PRIME_SAME'
				sameness_result_pos = the_5prime_edge_of_3prime_match

			# else: # sameness_result = 'DIFFERENT'

	return sameness_result, sameness_result_pos, sameness_result_pos2


######################################################
def record_unique_bam_sequence_strings():

	global ref_seq_chrom
	global ref_seq_start_pos
	global ref_seq_length
	global ref_seq_array
	global bam_read_chrom
	global bam_read_pos
	global bam_read_cigar
	global bam_read_sequence
	global bam_read_num_reads
	global bam_read_longest_read_length
	global seen_sequence
	global seen_sequence_count

	# Go through the bam reads twice.
	# For the first time, quickly look at only reads that match the ref.seq. to determine how much of ref.seq. has been seen, 
	# just by looking at the cigar string, without the more computationally expensive activity of looking at nucleotides, and save that seen ref.seq.
	# For the second time, look at only reads that don't match the ref.seq.
	# We need to look at and compare each nucleotide, which is computationally expensive.
	# Thus we avoid performing computationally expensive activities on reads that match the ref.seq. and thus don't need much inspection.

	# For the first time through looking at bam reads, quickly look at only reads that match the ref.seq. to determine how much of ref.seq. has been seen.

	seen_ref_seq_start_pos = []
	seen_ref_seq_end_pos = []
	for bam_upto in range( 0, bam_read_num_reads ):

		this_bam_cigar = bam_read_cigar[bam_upto]

		# is it a pattern such as 150M
		regex3 = re.compile(r'[0-9]+M$')
		regex_result = re.match( regex3, this_bam_cigar )
		if (regex_result):
			this_bam_start_pos = bam_read_pos[bam_upto]
			this_bam_end_pos = this_bam_start_pos + len(bam_read_sequence[bam_upto]) - 1
			this_seen_seq_is_recorded = False
			for i in range( 0, len(seen_ref_seq_start_pos) ):
				if ((this_bam_start_pos >= seen_ref_seq_start_pos[i]) and (this_bam_start_pos <= seen_ref_seq_end_pos[i])):
					seen_ref_seq_end_pos[i] = this_bam_end_pos
					this_seen_seq_is_recorded = True
				elif ((this_bam_end_pos >= seen_ref_seq_start_pos[i]) and (this_bam_end_pos <= seen_ref_seq_end_pos[i])):
					seen_ref_seq_start_pos[i] = this_bam_start_pos
					this_seen_seq_is_recorded = True
				elif ((this_bam_start_pos <= seen_ref_seq_start_pos[i]) and (this_bam_end_pos >= seen_ref_seq_end_pos[i])):
					seen_ref_seq_start_pos[i] = this_bam_start_pos
					seen_ref_seq_end_pos[i] = this_bam_end_pos
					this_seen_seq_is_recorded = True
			if (this_seen_seq_is_recorded == False):
				seen_ref_seq_start_pos.append( this_bam_start_pos )
				seen_ref_seq_end_pos.append( this_bam_end_pos )
		# record the ref_seq that we've seen in seen_sequence
		for this_seen_ref_seq_info in range( 0, len(seen_ref_seq_start_pos) ):
			this_bam_start_pos = seen_ref_seq_start_pos[this_seen_ref_seq_info]
			this_bam_end_pos = seen_ref_seq_end_pos[this_seen_ref_seq_info]
			this_bam_seq_aligned = [' '] * ref_seq_length
			this_bam_seq_aligned_count = [0] * ref_seq_length
			fill_start = this_bam_start_pos - ref_seq_start_pos
			fill_end_plus_1 = this_bam_end_pos - ref_seq_start_pos + 1
			ref_idx = 0
			if (fill_start >= 0):
				if ((fill_end_plus_1) > ref_seq_length):
					fill_end_plus_1 = ref_seq_length
				for i in range( fill_start, fill_end_plus_1 ):
					this_bam_seq_aligned[i] = ref_seq_array[ ref_idx : ref_idx+1 ]
					# Actually, we have probably seen this ref.seq. nucleotide more than once in the bam file.
					# We simply record as count = 1 because we don't want to spend computational time on ref.seq. bam reads
					# and read counts are important in our processing for non ref.seq. bam reads.
					this_bam_seq_aligned_count[i] = 1
					ref_idx = ref_idx + 1

	# For the second time through looking at bam reads, look at only reads that don't match the ref.seq.

	for bam_upto in range( 0, bam_read_num_reads ):

		this_bam_cigar = bam_read_cigar[bam_upto]
		# is it a pattern such as 150M
		regex4 = re.compile(r'[0-9]+M$')
		regex_result = re.match( regex4, this_bam_cigar )
		if (regex_result):
			do_nothing = True
		else:

			this_bam_seq = bam_read_sequence[bam_upto]
			this_bam_cigar = bam_read_cigar[bam_upto]
			this_bam_start_pos = bam_read_pos[bam_upto]

			this_bam_read_length = len(this_bam_seq)
			if (this_bam_read_length > bam_read_longest_read_length):
				bam_read_longest_read_length = this_bam_read_length

			this_bam_seq_aligned = [' '] * ref_seq_length
			this_bam_seq_aligned_count = [0] * ref_seq_length
			fill_start = this_bam_start_pos - ref_seq_start_pos
			fill_end_plus_1 = fill_start + len(this_bam_seq)
			bam_idx = 0
			if (fill_start < 0):
				bam_idx = fill_start * -1
				fill_start = 0
			if ((fill_end_plus_1) > ref_seq_length):
				fill_end_plus_1 = ref_seq_length
			for i in range( fill_start, fill_end_plus_1 ):
				this_bam_seq_aligned[i] = this_bam_seq[ bam_idx ]
				this_bam_seq_aligned_count[i] = 1
				bam_idx = bam_idx + 1

			save_5_PRIME_SAME_idx = -1
			save_5_PRIME_SAME_the_3prime_edge_of_5prime_match = -1
			save_3_PRIME_SAME_idx = -1
			save_3_PRIME_SAME_the_5prime_edge_of_3prime_match = -1
			save_MATCHING_EDGES_NOT_MIDDLE_idx = -1
			save_MATCHING_EDGES_NOT_MIDDLE_the_3prime_edge_of_5prime_match = -1
			save_MATCHING_EDGES_NOT_MIDDLE_the_5prime_edge_of_3prime_match = -1
			save_MATCHING_MIDDLE_idx = -1
			save_MATCHING_MIDDLE_the_3prime_edge_of_5prime_match = -1
			save_MATCHING_MIDDLE_the_5prime_edge_of_3prime_match = -1

			# We look at a given region of the chromosome.
			# If this bam read is barely in that region, then we would not be looking at much of it, so ignore it. Don't process only a bit of a read.

			if ( len(collapse_and_strip_seq(this_bam_seq_aligned)) > 50):

				# Look at each sequence.
				# Is it a continuation of a sequence we've already seen?
				# Or partway through the sequence, has it forked to become a new sequence?

				if (len(seen_sequence) == 0):

					# Add first bam sequence we see to the seen-concatenated-sequences list

					seen_sequence.append( this_bam_seq_aligned )
					seen_sequence_count.append( this_bam_seq_aligned_count )

				else:

					# See if this bam sequence overlap-matches sequences we've already seen

					already_seen_this_seq = False
					for seen_upto in range( 0, len(seen_sequence) ):
						if (already_seen_this_seq == False):

							this_seen_seq = seen_sequence[seen_upto]
							this_seen_seq_count = seen_sequence_count[seen_upto]

							sameness_result, sameness_result_pos, sameness_result_pos2 = are_these_2_sequences_the_same( this_bam_seq_aligned, this_seen_seq )

							if (sameness_result == 'SAME'):

								# This bam sequence is similar to a previously seen and stored concatenated-sequence, 
								# so update the seen-counts for that sequence.
								# If this bam sequence starts before any of the previously seen matching sequences, then fill in the start nucleotides.
								# Then write back this new state of the seen-concatenated-sequence to the storage array.

								already_seen_this_seq = True
								for i in range( 0, len(this_bam_seq_aligned) ):
									if (this_seen_seq[i] == ' '):
										this_seen_seq[i] = this_bam_seq_aligned[i]
									seen_sequence[seen_upto][i] = this_seen_seq[i]
									if (this_bam_seq_aligned[i] != ' '):
										if ((this_bam_seq_aligned[i] != 'N') or (seen_sequence[seen_upto][i] == ' ')):
											seen_sequence_count[seen_upto][i] = seen_sequence_count[seen_upto][i] + 1

							if (sameness_result == 'SAME_BUT_WITH_STUTTER'):

								# If we remove poly-A or poly-T sequences,
								# then this bam sequence is similar to a previously seen and stored concatenated-sequence, 
								# so update the seen-counts for that sequence.
								# Don't copy across any nucleotide sequence though. 
								# We don't know where the sequence should go. All we know is that it is stuttered and does notn't match sequence positions.
								# Do add to the counts to show that this part of the sequence has yet another read supporting its presence.

								already_seen_this_seq = True
								for i in range( 0, len(this_bam_seq_aligned) ):
									if (this_bam_seq_aligned[i] != ' '):
										seen_sequence_count[seen_upto][i] = seen_sequence_count[seen_upto][i] + 1

							elif (sameness_result == 'MATCHING_MIDDLE'):

								if (save_MATCHING_MIDDLE_idx == -1):
									save_MATCHING_MIDDLE_idx = seen_upto
									save_MATCHING_MIDDLE_the_3prime_edge_of_5prime_match = sameness_result_pos
									save_MATCHING_MIDDLE_the_5prime_edge_of_3prime_match = sameness_result_pos2

							elif (sameness_result == '5_PRIME_SAME'):

								if (save_5_PRIME_SAME_idx == -1):
									save_5_PRIME_SAME_idx = seen_upto
									save_5_PRIME_SAME_the_3prime_edge_of_5prime_match = sameness_result_pos

							elif (sameness_result == '3_PRIME_SAME'):

								if (save_3_PRIME_SAME_idx == -1):
									save_3_PRIME_SAME_idx = seen_upto
									save_3_PRIME_SAME_the_5prime_edge_of_3prime_match = sameness_result_pos

							elif (sameness_result == 'MATCHING_EDGES_NOT_MIDDLE'):

								if (save_MATCHING_EDGES_NOT_MIDDLE_idx == -1):
									save_MATCHING_EDGES_NOT_MIDDLE_idx = seen_upto
									save_MATCHING_EDGES_NOT_MIDDLE_the_3prime_edge_of_5prime_match = sameness_result_pos
									save_MATCHING_EDGES_NOT_MIDDLE_the_5prime_edge_of_3prime_match = sameness_result_pos2

					# If was saw a saved seen_sequence that was the 'SAME' as this bam_read_sequence,
					# then we would have added it to the list of saved seen_sequences already, we would have done it immediately.
					# However, if we saw any 5_PRIME_SAME, 3_PRIME_SAME, MATCHING_EDGES_NOT_MIDDLE, MATCHING_MIDDLE
					# we would have kept on looking at saved seen_sequences in the hopes of seeing one that is the SAME.
					# Now that we have finished comparing this bam_read_sequence to the saved seen_sequences,
					# if we didn't see any SAME, then look at any 5_PRIME_SAME, 3_PRIME_SAME, MATCHING_EDGES_NOT_MIDDLE, MATCHING_MIDDLE that we saw,
					# and add the new sequence as one of those (so that it gets the entire saved seen_sequence, not just the 150 bp of this bam_read_sequence).

					if (already_seen_this_seq == False):

						if (save_MATCHING_MIDDLE_idx != -1):

							seen_sequence.append( this_bam_seq_aligned )
							seen_sequence_count.append( this_bam_seq_aligned_count )

						elif (save_5_PRIME_SAME_idx != -1):

#							new_seq = []
							this_seen_seq = seen_sequence[save_5_PRIME_SAME_idx]
#							for i in range( 0, len(this_seen_seq) ):
#								new_seq.append( this_seen_seq[i] )
							new_seq = this_seen_seq
							the_3prime_edge_of_5prime_match = save_5_PRIME_SAME_the_3prime_edge_of_5prime_match
							for i in range( the_3prime_edge_of_5prime_match, len(this_bam_seq_aligned) ):
								new_seq[i] = this_bam_seq_aligned[i]
							new_seq_count = [0] * len(new_seq)
							for i in range( 0, len(new_seq) ):
								if (new_seq[i] != ' '):
									new_seq_count[i] = 1
							seen_sequence.append( new_seq )
							seen_sequence_count.append( new_seq_count )

						elif (save_3_PRIME_SAME_idx != -1):

							new_seq = []
							this_seen_seq = seen_sequence[save_3_PRIME_SAME_idx]
							for i in range( 0, len(this_seen_seq) ):
								new_seq.append( this_seen_seq[i] )
							the_5prime_edge_of_3prime_match = save_3_PRIME_SAME_the_5prime_edge_of_3prime_match
							for i in range( 0, (the_5prime_edge_of_3prime_match+1) ):
								new_seq[i] = this_bam_seq_aligned[i]
							new_seq_count = [0] * len(new_seq)
							for i in range( 0, len(new_seq) ):
								if (new_seq[i] != ' '):
									new_seq_count[i] = 1
							seen_sequence.append( new_seq )
							seen_sequence_count.append( new_seq_count )

						elif (save_MATCHING_EDGES_NOT_MIDDLE_idx != -1):

							seen_sequence.append( this_bam_seq_aligned )
							seen_sequence_count.append( this_bam_seq_aligned_count )

						elif (sameness_result == 'DIFFERENT'):

							seen_sequence.append( this_bam_seq_aligned )
							seen_sequence_count.append( this_bam_seq_aligned_count )

						else:
							raise ValueError("\n\nin Mobelwrap_identify_split_read_clusters.py, program logic error. Unexpected sameness_result. Thus will not continue processing any more sequences.\n")

	return

######################################################
def choose_final_list_of_unique_bam_sequence_strings( indexed_reference_genome_fasta_file, ref_seq_chrom, ref_seq_start_pos, ref_seq_length ):

	# Only keep concatenated bam sequences that we have seen more than once.
	# If the part of the sequence that is different to the reference genome has only been seen once,
	# then discard this sequence. 
	# If we haven't received a reference genome on input, then can't check that it is different to reference genome.

	global seen_sequence
	global seen_sequence_count
	global final_seen_sequence
	global final_seen_sequence_count
	global ref_seq_array

	minimum_number_of_times_seen = 1
	minimum_number_of_bp_to_have_been_seen_minimum_number_of_times = 10

	# look at each unique concatenated-bam-reads sequence to decide whether to keep it (because this sequence was seen more than once in the bam file) 
	# or not (because this sequence is from a read that was only seen once in the bam file)

	for seq_upto in range( 0, len(seen_sequence) ):

		this_seen_seq = seen_sequence[seq_upto]
		this_seen_seq_count = seen_sequence_count[seq_upto]
		keep_this_seq = False

		# look at the counts of how often the base pairs of this concatenated-bam-reads sequence were seen, in a sequence, in the bam file

		num_bp_seen_min_num_times = 0
		for i in range( 0, len(this_seen_seq_count) ):
			if (this_seen_seq_count[i] >= minimum_number_of_times_seen):
				num_bp_seen_min_num_times = num_bp_seen_min_num_times + 1
		if (num_bp_seen_min_num_times >= minimum_number_of_bp_to_have_been_seen_minimum_number_of_times):
			keep_this_seq = True

		# Compare this sequence to the reference genome.
		# If the novel parts of this sequence that do not match the genome have only been seen once,
		# then don't keep the sequence.
		# The rest of the sequence that matches the genome has probably been seen more than once.

		if ((keep_this_seq) and (indexed_reference_genome_fasta_file != '-')):

			accumulated_counter_for_novel_base_pairs = 0
			number_of_novel_base_pairs = 0
			number_of_base_pairs = 0
			length_to_compare = min( len(ref_seq_array), len(this_seen_seq) )

			for i in range( 0, length_to_compare ):
				if ((ref_seq_array[i] != 'N') and (this_seen_seq[i] != ' ')):
					number_of_base_pairs = number_of_base_pairs + 1
					if (ref_seq_array[i] != this_seen_seq[i]):
						accumulated_counter_for_novel_base_pairs = accumulated_counter_for_novel_base_pairs + this_seen_seq_count[i]
						number_of_novel_base_pairs = number_of_novel_base_pairs + 1
			similarity_of_this_sequence_to_reference_sequence = 0
			if (number_of_base_pairs > 0):
				similarity_of_this_sequence_to_reference_sequence = float(number_of_base_pairs - number_of_novel_base_pairs) / float(number_of_base_pairs)
				# if this concatenated-bam sequence is basically the same as the reference sequence,
				# then don't discard it when the few difference base pairs have only been seen once.
				if (similarity_of_this_sequence_to_reference_sequence < 0.9):
					average_counter_of_novel_base_pairs = float(accumulated_counter_for_novel_base_pairs) / float(number_of_novel_base_pairs)
					# if the novel parts of this concatenated-bam sequence have only been seen once,
					# then don't keep this concatenated-bam sequence.
					if (average_counter_of_novel_base_pairs <= 0.9):
						keep_this_seq = False

		if (keep_this_seq):
			final_seen_sequence.append( this_seen_seq )
			final_seen_sequence_count.append( this_seen_seq_count )

	return

######################################################
def main():

	global ref_seq_chrom
	global ref_seq_start_pos
	global ref_seq_length
	global bam_read_chrom
	global bam_read_pos
	global bam_read_cigar
	global bam_read_sequence
	global bam_read_num_reads
	global bam_read_longest_read_length
	global seen_sequence
	global seen_sequence_count
	global final_seen_sequence
	global final_seen_sequence_count

	ref_seq_chrom = ''
	ref_seq_start_pos = ''
	ref_seq_length = 0
	bam_read_chrom = []
	bam_read_pos = []
	bam_read_cigar = []
	bam_read_sequence = []
	bam_read_longest_read_length = 0
	seen_sequence = []
	seen_sequence_count = []
	final_seen_sequence = []
	final_seen_sequence_count = []

	### Read input parameters

	what_is_the_extract = sys.argv[1]
	where_is_bam = sys.argv[2]
	input_read_length = '-'
	read_length = 150
	if (len(sys.argv) >= 4):
		input_read_length = sys.argv[3]
		if (input_read_length != '-'):
			read_length = int(input_read_length)
	indexed_reference_genome_fasta_file = '-'
	if (len(sys.argv) >= 5):
		indexed_reference_genome_fasta_file = sys.argv[4]
	what_is_the_extract = what_is_the_extract.strip()
	where_is_bam = where_is_bam.strip()

	idx_colon = what_is_the_extract.rindex(':')
	idx_dash = what_is_the_extract.rindex('-')
	chrom = what_is_the_extract[0:idx_colon]
	chrom_first3chars = chrom[0:3]
	# chr_list = [ 'chr', 'Chr', 'CHR' ]
	# if (chrom_first3chars in chr_list):				# There is no need to remove the 'chr' because the input bam will have been aligned with the reference file, and both will either have 'chr' or not have it.
	# 	chrom = chrom[3:]
	chrom = str(chrom)
	input_num1 = int(what_is_the_extract[(idx_colon+1):idx_dash])
	input_num2 = int(what_is_the_extract[(idx_dash+1):])
	bam_num1 = input_num1 - read_length
	bam_num2 = input_num2 + read_length + read_length
	ref_num1 = bam_num1 - read_length
	ref_num2 = bam_num2 + read_length
	if (bam_num1 < 1):
		bam_num1 = 1
	if (bam_num2 < 1):
		bam_num2 = 1
	if (ref_num1 < 1):
		ref_num1 = 1
	if (ref_num2 < 1):
		ref_num2 = 1
	ref_seq_chrom = chrom
	ref_seq_start_pos = ref_num1
	ref_seq_length = ref_num2 - ref_num1 + 1

	# get the reference genome sequence for this part of the genome for which we are processing bam file reads

	if (indexed_reference_genome_fasta_file != '-'):
		command_output_lines = []
		ref_seq_end_pos = ref_seq_start_pos + ref_seq_length - 1
		samtools_faidx_position = ref_seq_chrom + ':' + str(ref_seq_start_pos) + '-' + str(ref_seq_end_pos)
		samtools_faidx_command = 'samtools faidx ' + indexed_reference_genome_fasta_file + ' ' + samtools_faidx_position
		command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
		if (command_status != 0):
			raise ValueError("\n\nIn Mobelwrap_identify_split_read_clusters.py, was not able to get the reference sequence from reference genome for this region of the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more sequences.\n")
		else:
			command_output_lines = parse_out_any_warnings( command_output, samtools_faidx_command )
		# some positions may not be found in reference genome, will have a header but no sequence
		if (len(command_output_lines) > 1):
			ref_seq_string = ''
			for i in range( 1, len(command_output_lines) ):
				ref_seq_string = ref_seq_string + command_output_lines[i]
			ref_seq_array = [' '] * len(ref_seq_string)
			for i in range( 0, len(ref_seq_string) ):
				this_bp = ref_seq_string[i]
				ref_seq_array[i] = this_bp

	# get all the reads in the bam file for this region

	get_bam_reads( where_is_bam, chrom, bam_num1, bam_num2 )
	if ((bam_read_longest_read_length > read_length) and (input_read_length == '-')):
		read_length = bam_read_longest_read_length
		bam_num1 = input_num1 - read_length
		bam_num2 = input_num2 + read_length
		if (bam_num1 < 1):
			bam_num1 = 1
		if (bam_num2 < 1):
			bam_num2 = 1
		if (ref_num1 < 1):
			ref_num1 = 1
		if (ref_num2 < 1):
			ref_num2 = 1
		ref_seq_chrom = chrom
		ref_seq_start_pos = ref_num1
		ref_seq_length = ref_num2 - ref_num1 + 1

	# record the strings of connected sequences seen in bam

	record_unique_bam_sequence_strings()

	# choose final list of strings of connected sequences seen in bam

	choose_final_list_of_unique_bam_sequence_strings( indexed_reference_genome_fasta_file, ref_seq_chrom, ref_seq_start_pos, ref_seq_length )

	# write the final list of unique concatenated bam sequences seen to STDOUT

	outline = str(ref_seq_chrom) + "\n"
	sys.stdout.write(outline)
	outline = str(ref_seq_start_pos) + "\n"
	sys.stdout.write(outline)
	for i in range( 0, len(final_seen_sequence) ):
		outline = collapse_seq( final_seen_sequence[i] ) + "\n"
		sys.stdout.write(outline)
		outline = collapse_seq_count( final_seen_sequence_count[i] ) + "\n"
		sys.stdout.write(outline)

if __name__=='__main__':
    main()


