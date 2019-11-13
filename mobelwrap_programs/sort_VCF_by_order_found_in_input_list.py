#!/usr/bin/python
# python sort_VCF_by_order_found_in_input_list.py -i input_VCF -l input_list -o output_VCF
# python sort_VCF_by_order_found_in_input_list.py -i MGRBp1_AABZR_sample1206_mobile_elements_withPolyACounts_sameSampleFields_sortedByChrom.vcf -l hs37d5x_sort_order.txt -o MGRBp1_AABZR_sample1206_mobile_elements_withPolyACounts_sameSampleFields_sortedByChrom_sortedByList.vcf

# This program fills the need to sort the VCF files by the same order as reference file hs37d5x.fa 
# so that GATK tools' CombineVariants can be used to merge multiple sample VCFs into 1 VCF having multiple samples.

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse

order_of_chromosomes = []
input_VCF_records = {}

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
	parser.add_argument('-l', action="store", dest="in_list", required=True, help='Input list containing one chromosome/contig per line, in the desired sort order')
	parser.add_argument('-o', action="store", dest="out_vcf", required=True, help='Output VCF file')
	args = parser.parse_args()

	# Read in the list containing the ordering of chromosomes/contigs.
	# Store the ordered list so that we can refer to it.

	in_list = open(args.in_list, 'r')
	for inline in in_list:
		inline = inline.strip()
		if (inline != ''):
			ref_chrom = str(inline)
			order_of_chromosomes.append( ref_chrom )

	# Read in the VCF file records.
	# Store in memory as an array (thus same as input order) within each chromosome,
	# with each new chromosome encountered in a dictionary identified by the chromosome,
	# to be written out later in a different order to the input order.
	# The VCF header is output immediately.

	out_vcf = open(args.out_vcf, 'w')
	in_vcf = open(args.in_vcf, 'r')
	current_chrom_vcf_records = []
	current_chrom = ''
	for inline in in_vcf:
		stripped_inline = inline.strip()
		if (stripped_inline != ''):
			first_char = inline[0:1]
			if (first_char == '#'):
				out_vcf.write( inline )
			else:
				infields = inline.split("\t")
				this_chrom = infields[0]
				if (this_chrom != current_chrom):
					# store the VCF records of the previous chromosome
					if (len(current_chrom_vcf_records) > 0):
						input_VCF_records[current_chrom] = current_chrom_vcf_records
					# start a new array for VCF records of this new chromosome
					current_chrom_vcf_records = []
					current_chrom = this_chrom
					# make sure this new chromosome is in the reference sort list
					# if it isn't, then add it to the ned
					found_this_chrom = False
					for ref_chrom in order_of_chromosomes:
						if (ref_chrom == this_chrom):
							found_this_chrom = True
					if (found_this_chrom == False):
						order_of_chromosomes.append( this_chrom )
				current_chrom_vcf_records.append( inline )
	if (len(current_chrom_vcf_records) > 0):
		input_VCF_records[current_chrom] = current_chrom_vcf_records

	# Write out the stored VCF records in the same chromosome order as the input list.

	for this_chrom in order_of_chromosomes:
		if this_chrom in input_VCF_records:
			this_chrom_vcf_records = input_VCF_records[this_chrom]
			for this_vcf_record in this_chrom_vcf_records:
				out_vcf.write( this_vcf_record )
	out_vcf.close()

if __name__=='__main__':
    main()


