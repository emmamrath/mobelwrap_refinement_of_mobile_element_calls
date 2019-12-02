#!/usr/bin/python
# python reformat_sam_file_aligned_against_ref.py input_sam input_fasta_ref
# python reformat_sam_file_aligned_against_ref.py CHM1_Illumina_all3_withM_sorted_markedDupl_noH_extract_1_104932100_104932250.sam ../GRCh37_fasta_extracts/GRCh37_1_104931000_104934000.fasta

# echo ">dna:chromosome chromosome:GRCh37:17:7652520:7652720" > temp.txt
# python reformat_sam_file_aligned_against_ref.py CHM1_Illumina_all3_withM_sorted_markedDupl_noH_extract_17_7652520_7652720.sam temp.txt


# The input sam file contains DNA read fragments, one per line, mapped to a genome.
# The input fasta file contains the part of the reference genome that the sam reads map to.
# The program outputs a file containing one line for the reference, 
# and one line for each sam read, lined up under the reference where it is mapped.
# The output file name is the input file name plus '.aligned.txt'

# Example input_sam:
# SRR1514952.182769727	83	1	104932009	60	101M	=	104931801	-309	ACAAAGACTTTTTAGATAAAAATCTCAAACATAAAAGCAGCAAAAGCAAAAATAGACAAAGGATTTCATCAAATTAAAAAGCTTCTGCACAAAAAAAGAAA	DDDEDDDDDCA@>CDDBEEDDDCFDDFFHGHHFHGIGCIIJJJJIGFJJJJIJJIJJJIJIHFBIICIEJJJIGJIJJJIIIGIJIIIHHGHGFFFFFCCC	MD:Z:1T99	PG:Z:MarkDuplicates.5	NM:i:1	AS:i:99	XS:i:26
# SRR1514950.139920169	99	1	104932011	60	101M	=	104932177	267	AAAGACTTTTTAGATAAAAATCTCAAACATAAAAGCAGCAAAAGCAAAAATAGACAAAGGATTTCATCAAATTAAAAAGCTTCTGCACAAAAAAAGAAAAC	@CCFFFFFHHHHHGIJJJJJIJIJIHGIJIIJIGCBDEEGIJIGAGGIIIJIGIGIJIJGHIEHGECEHIGEEGEHEFADC>BCEAECEEDDDDDDDDDD9	MD:Z:101	PG:Z:MarkDuplicates	NM:i:0	AS:i:101	XS:i:26
# SRR1514950.114999135	99	1	104932017	60	101M	=	104932191	275	TTTTTAGATAAAAATCTCAAACATAAAAGCAGCAAAAGCAAAAATAGACAAAGGATTTCATCAAATTAAAAAGCTTCTGCACAAAAAAAGAAAACAATCAA	CCCFFFFFHHHHHJJJJJJJJJJJJJJJIIJIFIJJIHJJJJJIIIJJJJJJJHIIGIJJJJJJJIJJJJIJIFHHHFFFFFFFEDCDDDDDDDDDDCCCC	MD:Z:101	PG:Z:MarkDuplicates	NM:i:0	AS:i:101	XS:i:26
# ST-E00197:144:HWLG3CCXX:1:2214:8603:471931   566263      16S109M25S  21  23080360

# Example input_fasta_ref:
# >1:9925-10525
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA
# ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA
# ACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTA


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
def main():

	### Read input files

	input_sam_file_name = sys.argv[1]
	input_fasta_ref_file_name = sys.argv[2]
	input_sam_file_name = input_sam_file_name.strip()
	input_fasta_ref_file_name = input_fasta_ref_file_name.strip()

	input_fasta_ref = open(input_fasta_ref_file_name, 'r')
	input_fasta_ref_lines = input_fasta_ref.readlines()
	input_fasta_ref.close()

	fasta_hdr = input_fasta_ref_lines[0]
	fasta_hdr = fasta_hdr[1:]
	fasta_hdr = fasta_hdr.replace('-',':')
	fasta_hdr_fields = fasta_hdr.split(':')
	ref_name = 'human'
	ref_chrom = fasta_hdr_fields[0]
	ref_start = int(fasta_hdr_fields[1])
	ref_end = int(fasta_hdr_fields[2])
	ref_sequence = ''
	for i in range( 1, len(input_fasta_ref_lines) ):
		inline = input_fasta_ref_lines[i]
		inline = inline.strip()
		ref_sequence = ref_sequence + inline

	input_sam = open(input_sam_file_name, 'r')
	input_sam_lines = input_sam.readlines()
	input_sam.close()

	outfile_name = input_sam_file_name + '.aligned.txt'
	outfile = open(outfile_name, 'w')

	### output the reference sequence

	outline_array = []
	reference_array = []

	ref_start_display = str(ref_start) + ':'
	leftcol = "%130s" % str(ref_start_display)
	grid = ''
	for i in range( ref_start, (ref_end + 1) ):
		grid_char = ' '
		if (((i % 10) == 7) or ((i % 10) == 8) or ((i % 10) == 9)):
			grid_char = ''
		if ((i % 10) == 0):
			grid_char = str(i)
			chr_len = len(grid_char)
			grid_char = grid_char[(chr_len-4):chr_len]
			if (i <= (ref_start+2)):
				chr_len = len(grid_char)
				if (i == ref_start):
					grid_char = grid_char[(chr_len-1):chr_len]
				elif (i == (ref_start+1)):
					grid_char = grid_char[(chr_len-2):chr_len]
				elif (i == (ref_start+2)):
					grid_char = grid_char[(chr_len-3):chr_len]
		grid = grid + grid_char
	outline = leftcol + grid + ':' + str(ref_end) + ' ' + "\n"
	reference_array.append( outline )

	leftcol = "%-130s" % ref_name
	grid = ''
	for i in range( ref_start, (ref_end + 1) ):
		grid_char = '-'
		if ((i % 5) == 0):
			grid_char = '+'
		if ((i % 10) == 0):
			grid_char = str(i)
			chr_len = len(grid_char)
			grid_char = grid_char[(chr_len-2):(chr_len-1)]
		grid = grid + grid_char
	outline = leftcol + grid + "\n"
	reference_array.append( outline )

	leftcol = "%-130s" % ref_name
	outline = leftcol + grid + "\n"
	outline = leftcol + ref_sequence + "\n"
	reference_array.append( outline )

	for refline in reference_array:
		outline_array.append( refline )

	### output the sam read alignments

	seq_line_count = 0
	for inline in input_sam_lines:
		inline = inline.strip()
		if (inline[0:1] != '@'):
			infields = inline.split("\t")
			qname = infields[0]
			flag = infields[1]
			rname = infields[2]
			pos = infields[3]
			mapq = infields[4]
			cigar = infields[5]
			rnext = infields[6]
			pnext = infields[7]
			tlen = infields[8]
			seq = infields[9]
			qual = infields[10]
			qname_pr = "%-54s" % qname
			rname_pr = "%-4s" % rname
			pos_pr = "%-12s" % pos
			cigar_pr = "%-12s" % cigar
			if (len(cigar_pr) > 12):
				cigar_pr = cigar_pr[0:12]
			rnext_pr = "%-4s" % rnext
			pnext_pr = "%-12s" % pnext
			leftcol = qname_pr + rname_pr + pos_pr + cigar_pr + rnext_pr + pnext_pr
			leftcol = leftcol
			leftcol = "%-130s" % leftcol
			first_soft_clip = 0
			cigar_split = cigar.split('S')
			if (is_integer(cigar_split[0])):
				first_soft_clip = int(cigar_split[0])
			padding_needed = int(pos) - ref_start - first_soft_clip
			if (padding_needed < 0):
				padding_needed = 0
			padding = ' ' * padding_needed
			outline = leftcol + padding + seq + "\n"
			outline_array.append( outline )

			# output the reference sequence every 50 lines of sam fragments
			seq_line_count = seq_line_count + 1
			if ((seq_line_count % 50) == 0):
				for refline in reference_array:
					outline_array.append( refline )
				seq_line_count = 0

	for outline in outline_array:
		# sys.stdout.write(outline)
		outfile.write(outline)

	outfile.write("\n\n")
	outfile.close()


if __name__=='__main__':
    main()

