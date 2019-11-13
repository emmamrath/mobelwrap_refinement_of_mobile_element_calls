#!/usr/bin/python
# python prepare_VCF_after_running_GATK.py -i input_VCF -o output_VCF
# python prepare_VCF_after_running_GATK.py -i MGRBp1_AABZR_sample1206_mobile_elements_withPolyACounts_sameSampleFields_sortedByChrom_sortedByList.vcf -o MGRBp1_AABZR_sample1206_mobile_elements_withPolyACounts_sameSampleFields_sortedByChrom_sortedByList_pos1.vcf

# This program changes back various fields in a VCF file after it was run using GATK tools.
#   *   For any MEI format fields in the sample fields that are set to a dash to signify no value, change it back to a dot.

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
def main():

	parser = argparse.ArgumentParser(description='Read in VCF file and input-list file. Sort the VCF by the chromosomes/contigs as they appear in the input list.')
	parser.add_argument('-i', action="store", dest="in_vcf", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="out_vcf", required=True, help='Output VCF file')
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

				# For any MEI format fields in the sample fields that are set to a dash to signify no value, change it back to a dot.

				if (len(infields) >= 10):
					for i in range( 9, len(infields) ):
						this_sample = infields[i]
						new_sample = ''
						this_sample_fields = this_sample.split(':')
						for this_sample_field in this_sample_fields:
							if (this_sample_field == '-'):
								this_sample_field = '.'
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


