#!/usr/bin/python
# python write_out_a_file_with_no_final_carriage_return.py input_file output_file

# python write_out_a_file_with_no_final_carriage_return.py NA12878_CHROM_13_redo1LinePerAllele_Cato_10variants.vcf NA12878_CHROM_13_redo1LinePerAllele_Cato_10variants_missingEnd.vcf

# This program creates a file that can be used for testing
# another program that needs to detect a file that has not been completely written out.

import sys
import os

input_file_name = sys.argv[1]
output_file_name = sys.argv[2]

infile = open( input_file_name, "r" )
inlines = infile.readlines()
infile.close()

outfile = open( output_file_name, "w" )

for i in range( 0, (len(inlines) - 1) ):
	inline = inlines[i]
	outfile.write( inline )

i = len(inlines) - 1
inline = inlines[i]
inline = inline.strip()
outfile.write( inline )

outfile.close()

