#!/usr/bin/python
# python clean_up_vcf_file_for_one_sample.py -i input_VCF -r reference_sequence -o output_VCF
# python clean_up_vcf_file_for_one_sample.py -i mei_AAAAA_crop.vcf -r /g/data3/zc9/resources/reference_genomes/hs37d5x/hs37d5x.fa -o mei_AAAAA_crop_cleaned.vcf
# python clean_up_vcf_file_for_one_sample.py -i bicseq_BAAWU_crop.vcf -r /g/data3/zc9/resources/reference_genomes/hs37d5x/hs37d5x.fa -o bicseq_BAAWU_crop_cleaned.vcf
# python clean_up_vcf_file_for_one_sample.py -i gridss_AACER_crop.vcf -r /g/data3/zc9/resources/reference_genomes/hs37d5x/hs37d5x.fa -o gridss_AACER_crop_cleaned.vcf
# python clean_up_vcf_file_for_one_sample.py -i gridss_AACER_crop.vcf -r /g/data3/zc9/resources/reference_genomes/hs37d5x/hs37d5x.fa -o gridss_AACER_crop_cleaned_gt.vcf -gt GT

# This program makes the following changes to a one-sample VCF file:
#   *   Set POS to 1 when POS is 0 or less than 1, and fill REF with the correct nucleotide for POS 1
#   *   If REF is a dot, then set it to the ref seq nucleotide for that POS, or to N if the chromosome or contig doesn't have that POS
#   *   If the sample GT field is a dot, set it to 1/.
#   *   If the REF is 1 bp and is not upper or lower case A,C,G,T,N, then set it to N. Otherwise GATK will crash if it sees an r or R (which was seen in an hs37d5x contig)


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse
import commands


######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def convert_to_integer(input_string):
	s_int = 0
	if ((input_string != 'NA') and (input_string != '')):
		s_int = int(input_string)
        return s_int

######################################################
def is_valid_ref(input_ref):
	is_valid = False
	valid_chars = 'ACGTNacgtn'
	if (valid_chars.find(input_ref) >= 0):
		is_valid = True
        return is_valid

######################################################
def main():

	parser = argparse.ArgumentParser(description='Read in VCF file, fix a few fields, output the resulting VCF file.')
	parser.add_argument('-i', action="store", dest="in_vcf", required=True, help='Input VCF file')
	parser.add_argument('-r', action="store", dest="in_ref_seq", required=True, help='Input reference sequence in fasta format and indexed so it can be called by samtools')
	parser.add_argument('-o', action="store", dest="out_vcf", required=True, help='Output VCF file')
	parser.add_argument('-gt', action="store", dest="gt_flag", required=False, help='If this gt flag is present, do not change the GT field. Otherwise, GT field of dot will be set to 1/.')
	args = parser.parse_args()

	# Read in and process each VCF file record.
	# Then write out the record.

	out_vcf = open(args.out_vcf, 'w')
	in_vcf = open(args.in_vcf, 'r')
	for inline in in_vcf:
		stripped_inline = inline.strip()
		if (stripped_inline != ''):
			first_char = inline[0:1]
			if (first_char == '#'):
				out_vcf.write( inline )
			else:
				infields = inline.split("\t")

				# Set POS to 1 when POS is 0 or less than 1, and fill REF with the correct nucleotide for POS 1
				# If REF is a dot, then set it to the ref seq nucleotide for that POS, or to N if the chromosome or contig doesn't have that POS

				this_chrom = infields[0]
				this_pos = int(infields[1])
				this_ref = infields[3]
				if (this_ref == '.'):
					infields[3] = 'N'
				if (len(this_ref) == 1):
					if (is_valid_ref(this_ref) == False):
						infields[3] = 'N'
				if ((this_pos < 1) or (this_ref == '.')):
					if (this_pos < 1):
						this_pos = 1
						infields[1] = 1
						infields[3] = 'N'

					samtools_faidx_position = this_chrom + ':' + str(this_pos) + '-' + str(this_pos)
					samtools_faidx_command = 'samtools faidx ' + args.in_ref_seq + ' ' + samtools_faidx_position
					command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
					if (command_status != 0):
						raise ValueError("\n\nWas not able to get the reference sequence from reference genome for this region of the genome using command:\n" + samtools_faidx_command + "\nThus will not continue processing any more VCF records.\n")
					command_output_lines = command_output.split("\n")
					# some positions may not be found in reference genome, will have a header but no sequence
					if (len(command_output_lines) > 1):
						infields[3] = command_output_lines[1].strip()

				# If GT field is a dot, then change it to 1/.

				this_format = infields[8]
				bits = this_format.split(':')
				this_format_1 = bits[0]
				if (this_format_1 == 'GT'):
					if (len(infields) >= 10):
						for i in range( 9, len(infields) ):
							this_sample = infields[i]
							new_sample = ''
							this_sample_fields = this_sample.split(':')
							for field_idx in range( 0, len(this_sample_fields) ):
								this_sample_field = this_sample_fields[field_idx]
								if (field_idx == 0): # this is the GT field
									if (this_sample_field == '.'):
										if (args.gt_flag is None):
											this_sample_field = '1/.'
										# else when the args.gt_flag is present, we don't touch the GT field, we leave it as . or ./.
								if (new_sample == ''):
									new_sample = str(this_sample_field)
								else:
									new_sample = new_sample + ':' + str(this_sample_field)
							infields[i] = new_sample

				outline = infields[0]
				for i in range( 1, len(infields) ):
					outline = outline + "\t" + str(infields[i])
				out_vcf.write( outline )
	out_vcf.close()

if __name__=='__main__':
    main()

