#!/usr/bin/python
# python set_qual_field_of_multi_sample_VCF.py 
# cat test1.vcf | python set_qual_field_of_multi_sample_VCF.py > test1_qual.vcf

# This program takes 1 VCF file on input from stdin, and produces 1 VCF file on output to stdout.

# This program reads the sample QUAL or MEIQUAL fields and set the VCF.QUAL record value to the highest of the sample values

# This adjustment to QUAL is necessary because GATK GenomeAnalysisTK.jar CombineVariants 
# fills VCF.QUAL with the VCF.QUAL of the first sample it met for that position, not the highest QUAL.

import sys
import os
import commands
import re

DELIMITER_FOR_INFO_MEINFO_FIELD = ','

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def is_float(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

######################################################
def choose_qual( qual1, qual2 ):

	qual1 = int(qual1)
	qual2 = int(qual2)
	return_qual = qual1
	if (qual2 > qual1):
		return_qual = qual2
        return return_qual

######################################################
def extract_sample_qual_index( this_inline_format ):

	return_qual_index = -1
	infields = this_inline_format.split(":")
	for i in range( 0, len(infields) ):
		this_field = infields[i]
		if (this_field == 'QUAL'):
			return_qual_index = i
		elif (this_field == 'MEIQUAL'):
			if (return_qual_index == -1): # QUAL takes precedence over MEIQUAL
				return_qual_index = i
        return return_qual_index

######################################################
def extract_qual_from_sample( qual_index, sample_fields_string ):

	return_qual = float(0)
	infields = sample_fields_string.split(":")
	return_qual_str = str(infields[qual_index])
	if (return_qual_str != '.'):
		return_qual = float(return_qual_str)
        return return_qual

######################################################
def main():

	# Read in the input VCF file from STDIN

	in_header = True
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == True):
			sys.stdout.write( inline )

		else: # We are processing VCF data records. We are no longer in the header part of the file.

			inline_stripped = inline.strip()
			infields = inline_stripped.split("\t")
			this_inline_chrom = str(infields[0])
			this_inline_pos = str(infields[1])
			this_inline_id = str(infields[2])
			this_inline_ref = str(infields[3])
			this_inline_alt = str(infields[4])
			this_output_qual = 0
			if (is_float(infields[5])):
				this_output_qual = float(infields[5])
			this_inline_filter = str(infields[6])
			this_inline_info = str(infields[7])
			this_inline_format = str(infields[8])
			sample_qual_index = extract_sample_qual_index( this_inline_format )

			all_output_samples = ''
			if (len(infields) > 9):
				for i in range( 9, len(infields) ):

					this_sample_inline = infields[i].strip()
					all_output_samples = all_output_samples + "\t" + this_sample_inline

					if (sample_qual_index > -1):
						if ((this_sample_inline != './.') and (this_sample_inline != '.')):
							this_sample_qual = extract_qual_from_sample( sample_qual_index, this_sample_inline )
							if (this_sample_qual > 0):
								this_output_qual = choose_qual( this_sample_qual, this_output_qual )

			outline = this_inline_chrom + "\t" + this_inline_pos + "\t" + this_inline_id + "\t" + this_inline_ref + "\t" + this_inline_alt
			outline = outline + "\t" + str(this_output_qual) + "\t" + this_inline_filter + "\t" + this_inline_info + "\t" + this_inline_format
			outline = outline + all_output_samples + "\n"
			sys.stdout.write( outline )


if __name__=='__main__':
    main()


