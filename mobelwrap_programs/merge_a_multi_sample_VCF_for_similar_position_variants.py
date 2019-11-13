#!/usr/bin/python
# python merge_a_multi_sample_VCF_for_similar_position_variants.py input_Mobster input_window_in_basepairs input_reference_genome_fastq
# cat MGRBp1_sample1000_sample630_merged_Mobster_predictions_sorted.vcf | python merge_a_multi_sample_VCF_for_similar_position_variants.py 100 /nvme/mobile_elements_2016dec/1000_genomes_hs37d5_fastq/human_g1k_v37.fasta > MGRBp1_sample1000_sample630_merged_Mobster_predictions_sorted_merged100bp.vcf

# This program takes 1 VCF file on input from stdin, and produces 1 VCF file on output to stdout.
# This program merges IMPRECISE variants having the same alleles 
# whose ALT allele starts with <INS:ME (thus they are mobile element insertion variants),
# whose positions are within input_window_in_basepairs of the first variant seen in the group of close variants.
# A chain of variants where the first and last variants are further aparts than input_window_in_basepairs will thus not be merged.
# The final POS of merged variants will be the average POS of the merged variants.
# This means that we have to get the reference sequence's nucleotide at that new POS position.

# Which field to take for QUAL?
# This program has the same behaviour as bcftools merge, that is,
# QUAL receives the highest QUAL value of the merged VCF records.
# Which fields to take for ID and FILTER?
# Let's take the fields from the first VCF record in the merge group (until we find a reason to do otherwise).
# Which fields to take for FORMAT?
# Let's take the field from the first VCF record, and verify that all subsequent records have the same FORMAT.
# If they don't, crash the program.

# Which fields to take for INFO?
# Within the INFO field, there is an MEINFO field with subfields START and END, and a TARGETSITEDUPL field.
# The START and END fields will be adjusted to contain the widest range of values seen on input.
# When the MEI variant has become multi-allele (eg. ALT=<INS:ME:L1>,<INS:ME:ALU> at the same POS),
# then MEINFO field for each allele is listed, comma-separated as for ALT field,
# and TARGETSITEDUPL field for each allele is listed, comma-separated.
# The other 2 fields in INFO are always IMPRECISE and INS, and they are taken from the first record seen (which will have values IMPRECISE and INS).

# The Mobster variants pipeline after Mobsters has called mobile element insertions from a bam file is:
# The pipeline uses vcftools: 	http://vcftools.sourceforge.net/perl_module.html   	git clone https://github.com/vcftools/vcftools.git
# The pipeline uses bcftools: 	https://samtools.github.io/bcftools/   			git clone git://github.com/samtools/bcftools.git bcftools
# The bcftools uses htslib: 	http://vcftools.sourceforge.net/htslib.html   		git clone --branch=develop git://github.com/samtools/htslib.git htslib
#
# python create_VCF_annotation_from_Mobster_predictions_for_one_sample.py Sample1_Mobster_predictions.txt human_reference_genome_indexed.fasta sample_mapping.tsv
# python create_VCF_annotation_from_Mobster_predictions_for_one_sample.py Sample2_Mobster_predictions.txt human_reference_genome_indexed.fasta sample_mapping.tsv
# python create_VCF_annotation_from_Mobster_predictions_for_one_sample.py Sample3_Mobster_predictions.txt human_reference_genome_indexed.fasta sample_mapping.tsv
# ...
# vcf-sort Sample1_Mobster_predictions.txt.vcf > Sample1_Mobster_predictions_sorted.vcf
# vcf-sort Sample2_Mobster_predictions.txt.vcf > Sample2_Mobster_predictions_sorted.vcf
# vcf-sort Sample3_Mobster_predictions.txt.vcf > Sample3_Mobster_predictions_sorted.vcf
# ...
# bgzip -c Sample1_Mobster_predictions_sorted.vcf > Sample1_Mobster_predictions_sorted.vcf.gz
# bgzip -c Sample2_Mobster_predictions_sorted.vcf > Sample2_Mobster_predictions_sorted.vcf.gz
# bgzip -c Sample3_Mobster_predictions_sorted.vcf > Sample3_Mobster_predictions_sorted.vcf.gz
# ...
# bcftools index Sample1_Mobster_predictions_sorted.vcf.gz
# bcftools index Sample2_Mobster_predictions_sorted.vcf.gz
# bcftools index Sample3_Mobster_predictions_sorted.vcf.gz
# ...
# bcftools merge -m all -o AllSamples_Mobster_predictions.vcf Sample1_Mobster_predictions_sorted.vcf.gz Sample2_Mobster_predictions_sorted.vcf.gz Sample3_Mobster_predictions_sorted.vcf.gz
# 
# cat AllSamples_Mobster_predictions.vcf | python merge_a_multi_sample_VCF_for_similar_position_variants.py 100 /nvme/mobile_elements_2016dec/1000_genomes_hs37d5_fastq/human_g1k_v37.fasta > AllSamples_Mobster_predictions_merged100bp.vcf


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
def decide_if_variant_is_type_can_be_merged( alt_field ):

	can_be_merged = False
	if (len(alt_field) >= 8):
		alt_substr = alt_field[0:8]
		if ((alt_substr == '<INS:ME:') or (alt_substr == '<DEL:ME')):
			can_be_merged = True
        return can_be_merged

######################################################
def decide_if_same_variant( window, chrom1, chrom2, pos1, pos2_array ):

	same_variant = False
	if (chrom1 == chrom2):
		pos1 = int(pos1)
		pos2 = ''
		if (len(pos2_array) > 0):
			pos2 = int(pos2_array[0])
		pos2 = int(pos2)
		if ( abs(pos1 - pos2) <= window ):
			same_variant = True
        return same_variant

######################################################
def choose_qual( qual1, qual2 ):

	qual1 = float(qual1)
	qual2 = float(qual2)
	return_qual = qual1
	if (qual2 > qual1):
		return_qual = qual2
        return return_qual

######################################################
def choose_best_meinfo( meinfo1, meinfo2 ): # eg. L1,-179,180,.

	global DELIMITER_FOR_INFO_MEINFO_FIELD

	new_meinfo = [ '', 0, 0, '' ] # MEI-type , start-pos , end-pos , polarity
	new_meinfo_string = ''
	if ((meinfo1 == '') and (meinfo2 == '')):
		new_meinfo_string = ''
	elif (meinfo1 == ''):
		new_meinfo_string = meinfo2
	elif (meinfo2 == ''):
		new_meinfo_string = meinfo1
	else:
		meinfo1_fields = meinfo1.split( DELIMITER_FOR_INFO_MEINFO_FIELD )
		meinfo2_fields = meinfo2.split( DELIMITER_FOR_INFO_MEINFO_FIELD )
		new_meinfo[0] = meinfo1_fields[0]
		new_meinfo[1] = int(meinfo1_fields[1])
		if (int(meinfo2_fields[1]) < int(meinfo1_fields[1])):
			new_meinfo[1] = int(meinfo2_fields[1])
		new_meinfo[2] = int(meinfo1_fields[2])
		if (int(meinfo2_fields[2]) > int(meinfo1_fields[2])):
			new_meinfo[2] = int(meinfo2_fields[2])
		new_meinfo[3] = meinfo1_fields[3]
		if ((meinfo1_fields[3] == '.') and (meinfo2_fields[3] != '.')):
			new_meinfo[3] = meinfo2_fields[3]
		new_meinfo_string = new_meinfo[0] + DELIMITER_FOR_INFO_MEINFO_FIELD + str(new_meinfo[1]) + DELIMITER_FOR_INFO_MEINFO_FIELD + str(new_meinfo[2]) + DELIMITER_FOR_INFO_MEINFO_FIELD + new_meinfo[3]
        return new_meinfo_string

######################################################
def choose_best_targetsitedupl( targetsitedupl1, targetsitedupl2 ): # values are duplication, noTSD, or unknown

	new_targetsitedupl = targetsitedupl1
	if (targetsitedupl1 == 'unknown'):
		new_targetsitedupl = targetsitedupl2
	elif ((targetsitedupl1 == 'noTSD') and (targetsitedupl2 == 'duplication')):
		new_targetsitedupl = 'unknown'
	elif ((targetsitedupl2 == 'noTSD') and (targetsitedupl1 == 'duplication')):
		new_targetsitedupl = 'unknown'
        return new_targetsitedupl

######################################################
def create_new_format_from_2_different_formats( format1, format2  ):

	new_format = ''
	got_GT = False
	format1_fields = format1.split(':')
	format2_fields = format2.split(':')
	new_format_fields = {}
	for this_field in format1_fields:
		new_format_fields[this_field] = True
	for this_field in format2_fields:
		if this_field not in new_format_fields:
			new_format_fields[this_field] = True
	for this_field in new_format_fields:
		if (this_field == 'GT'):
			got_GT = True
	if (got_GT):
		new_format = 'GT'
	for this_field in new_format_fields:
		if ((this_field == 'GT') and (got_GT)):
			do_nothing = True # We already added this GT genotype field as the first field in the new format
		else:
			if (new_format == ''):
				new_format = str(this_field)
			else:
				new_format = new_format + ':' + str(this_field)
        return new_format

######################################################
def convert_sample_to_new_format( one_sample, old_format, new_format ):

	# It would be nice to assume that all VCF records have the same FORMAT fields in the sample fields
	# because when merging 2 VCF records, they will get the same FORMAT field in the one record, and all samples need to have all those fields.
	# However, it appears that GATK CombineVariants removes some of the FORMAT fields, and its output is input to this program,
	# so we can't assume that all VCF records to be merged have the same FORMAT fields.

	one_sample_reformatted = ''
	if (one_sample == './.'):
		one_sample_reformatted = './.'
	else:
		one_sample_fields = one_sample.split(':')
		old_format_fields = old_format.split(':')
		new_format_fields = new_format.split(':')
		one_sample_reformatted_fields = {}
		for this_field_name in new_format_fields:
			one_sample_reformatted_fields[this_field_name] = '.'
		for i in range( 0, len(old_format_fields) ):
			this_field_name = old_format_fields[i]
			this_field_value = '.'
			if (len(one_sample_fields) > i):
				this_field_value = one_sample_fields[i]
			one_sample_reformatted_fields[this_field_name] = this_field_value
		for this_field_name in new_format_fields:
			this_field_value = one_sample_reformatted_fields[this_field_name]
			if (one_sample_reformatted == ''):
				one_sample_reformatted = str(this_field_value)
			else:
				one_sample_reformatted = one_sample_reformatted + ':' + str(this_field_value)

        return one_sample_reformatted

######################################################
def find_first_separator_position( this_field ):

	separator_position = -1
	separator = ''
	slash_position = this_field.find('/')
	bar_position = this_field.find('|')
	if (slash_position >= 0):
		separator_position = slash_position
		separator = '/'
	if (bar_position >= 0):
		if (bar_position < separator_position):
			separator_position = bar_position
			separator = '|'

	return separator_position, separator

######################################################
def main():

	input_window_in_basepairs = sys.argv[1]
	input_window_in_basepairs = int(input_window_in_basepairs)
	input_reference_genome_fastq_file = sys.argv[2] # needs to have been indexed with bwa index

	# Read in the input VCF file from STDIN

	in_header = True
	previous_ME_inline_chrom = ''
	previous_ME_inline_pos_array = []
	previous_ME_inline_id = ''
	previous_ME_inline_ref = ''
	previous_ME_inline_alt = ''
	previous_ME_inline_qual = ''
	previous_ME_inline_filter = ''
	previous_ME_inline_info = ''
	previous_ME_inline_format = ''
	previous_ME_inline_multiple_samples_fields = ''
	previous_ME_inline_subsequent_non_ME_inlines_to_write_out = []
	
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == True):
			sys.stdout.write( inline )

		else: # We are processing VCF data records. We are no longer in the header part of the file.

			# Is this variant of a type that we are interested in merging?

			inline_stripped = inline.strip()
			infields = inline_stripped.split("\t")
			this_inline_chrom = infields[0]
			this_inline_pos = infields[1]
			this_inline_id = infields[2]
			this_inline_ref = infields[3]
			this_inline_alt = infields[4]
			this_inline_qual = infields[5]
			this_inline_filter = infields[6]
			this_inline_info = infields[7]
			this_inline_format = infields[8]
			this_inline_multiple_samples_fields = ''
			if (len(infields) >= 10):
				this_inline_multiple_samples_fields = infields[9]
				if (len(infields) >= 11):
					for i in range( 10, len(infields) ):
						this_inline_multiple_samples_fields = this_inline_multiple_samples_fields + "\t" + infields[i]

			variant_is_type_that_can_be_merged_result = decide_if_variant_is_type_can_be_merged( this_inline_alt )

			if (variant_is_type_that_can_be_merged_result == False):

				# If this variant is not of a type that we are interested in merging,
				# then it will be written out without being merged.
				# It may need to written out after merging of a previous variant prior to it
				# which may not yet be ready and completely merged to write out,
				# so that the output file is still sorted.

				if (previous_ME_inline_chrom == ''):
					sys.stdout.write( inline )
				else:
					previous_ME_inline_subsequent_non_ME_inlines_to_write_out.append( inline )

			else: # This variant is of a type that we will merge if it's POS is close to similar variants.

				similar_position_result = decide_if_same_variant( input_window_in_basepairs, this_inline_chrom, previous_ME_inline_chrom, this_inline_pos, previous_ME_inline_pos_array )

				if (similar_position_result == False):

					# If this variant is not close enough to be merged with previous variant(s)
					# then write out any previous variants, 
					# and write out any non-ME variants being saved to write out after the previous merged variants.
					# Don't write out this variant though,
					# because it might be close to future variants and be merged with them.

					# Write out any previous ME variants being merged

					if (previous_ME_inline_chrom != ''):

						# When writing out a merged variant, get average POS and it's corresponding REF

						previous_ME_inline_avg_pos = previous_ME_inline_pos_array[0]

						if (len(previous_ME_inline_pos_array) > 1):

							# The previous variant is a bunch of variants of similar position.
							# Now that we are writing them out as one merged variant
							# and know all the variants in this merged bunch,
							# get the average POS and corresponding REF

							# Get the average POS

							sum_of_pos = 0
							for this_pos in previous_ME_inline_pos_array:
								sum_of_pos = sum_of_pos + this_pos
							previous_ME_inline_avg_pos = sum_of_pos / len(previous_ME_inline_pos_array)
							previous_ME_inline_avg_pos = int(previous_ME_inline_avg_pos)
							if (previous_ME_inline_avg_pos < 1):
								previous_ME_inline_avg_pos = 1

							# Get the corresponding REF of that POS

							samtools_faidx_position = str(previous_ME_inline_chrom) + ':' + str(previous_ME_inline_avg_pos) + '-' + str(previous_ME_inline_avg_pos)
							samtools_faidx_command = 'samtools faidx ' + input_reference_genome_fastq_file + ' ' + samtools_faidx_position
							command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
							if (command_status != 0):
								raise ValueError("\n\nWas not able to get REF position from reference genome. Thus will not continue processing any more mobile elements.\n")
							command_output_lines = command_output.split("\n")
							# some positions that mobile element mapped to may not be found in reference genome
							# eg. hs37d5:180728
							if (len(command_output_lines) > 1):
								previous_ME_inline_ref = command_output_lines[1].strip()

						outline = str(previous_ME_inline_chrom) + "\t" + str(previous_ME_inline_avg_pos) + "\t" + str(previous_ME_inline_id) + "\t" + previous_ME_inline_ref + "\t" + previous_ME_inline_alt + "\t" + str(previous_ME_inline_qual) + "\t" + previous_ME_inline_filter + "\t" + previous_ME_inline_info + "\t" + previous_ME_inline_format + "\t" + previous_ME_inline_multiple_samples_fields + "\n"
						sys.stdout.write( outline )

					# Write out any previous non-ME variants being saved to write out after the previous merged variants.

					if ( len(previous_ME_inline_subsequent_non_ME_inlines_to_write_out) > 0 ):
						for previous_ME_inline_subsequent_inline in previous_ME_inline_subsequent_non_ME_inlines_to_write_out:
							sys.stdout.write( previous_ME_inline_subsequent_inline )
						previous_ME_inline_subsequent_non_ME_inlines_to_write_out = []

					# Don't write out this ME variant though.
					# Save it so it can be merged with future ME variants if need be.

					previous_ME_inline_chrom = this_inline_chrom
					previous_ME_inline_pos_array = []
					previous_ME_inline_pos_array.append( int(this_inline_pos) )
					previous_ME_inline_id = this_inline_id
					previous_ME_inline_ref = this_inline_ref
					previous_ME_inline_alt = this_inline_alt
					previous_ME_inline_qual = this_inline_qual
					previous_ME_inline_filter = this_inline_filter
					previous_ME_inline_info = this_inline_info
					previous_ME_inline_format = this_inline_format
					previous_ME_inline_multiple_samples_fields = this_inline_multiple_samples_fields

				else: 

					# This ME variant needs to be merged with previous ME variant(s). 
					# There is at least one previous ME variant to merge this one with
					# because we have just seen that this ME variant is within the similarity window
					# of the first of the previous ME variants in the currently saved previous ME variants bunch.

					# Verify that variants to be merged have the same FORMAT field.
					# If they don't, then make them the same by expanding both to have all fields from both.

					if (previous_ME_inline_format != this_inline_format):
						new_inline_format = create_new_format_from_2_different_formats( this_inline_format, previous_ME_inline_format  )
						new_samples_fields = ''
						this_inline_multiple_samples_fields_split = this_inline_multiple_samples_fields.split("\t")
						i = 1
						for one_sample in this_inline_multiple_samples_fields_split:
							one_sample_reformatted = convert_sample_to_new_format( one_sample, this_inline_format, new_inline_format )
							if (new_samples_fields == ''):
								new_samples_fields = one_sample_reformatted
							else:
								new_samples_fields = new_samples_fields + "\t" + one_sample_reformatted
							i = i + 1
						this_inline_multiple_samples_fields = new_samples_fields
						new_samples_fields = ''
						previous_ME_inline_multiple_samples_fields_split = previous_ME_inline_multiple_samples_fields.split("\t")
						for one_sample in previous_ME_inline_multiple_samples_fields_split:
							one_sample_reformatted = convert_sample_to_new_format( one_sample, previous_ME_inline_format, new_inline_format )
							if (new_samples_fields == ''):
								new_samples_fields = one_sample_reformatted
							else:
								new_samples_fields = new_samples_fields + "\t" + one_sample_reformatted
						previous_ME_inline_multiple_samples_fields = new_samples_fields
						this_inline_format = new_inline_format
						previous_ME_inline_format = new_inline_format

					# Choose the highest QUAL from amongst the variants to be merged

					previous_ME_inline_qual = choose_qual( previous_ME_inline_qual, this_inline_qual )

					# Keep a list of all the POS being merged. Will choose the average when we have all the POS to be merged.

					previous_ME_inline_pos_array.append( int(this_inline_pos) )

					# If the new VCF record's ALTs are different to the existing ALTs,
					# then add the new ALTs to the existing ALTS.
					# Then add each individual new sample's formatted fields 
					# as new sample-formatted-field columns
					# such that any Genotype GT field values point to the new ALT value and not to their old ALT value. 

					# Make the new ALT field as a string list of all the ALTs in the the previous ME variant
					# plus any new unique ALT values in this new VCF record.

					# merging of MEINFO and TARGETSITEDUPL in INFO field, eg.
					# if merging 2 ALT to become <INS:ME:ALU>,<INS:ME:L1> then their INFO fields:
					# IMPRECISE;SVTYPE=INS;MEINFO=L1,-179,180,.;TARGETSITEDUPL=unknown
					# IMPRECISE;SVTYPE=INS;MEINFO=ALU,-16,16,.;TARGETSITEDUPL=duplication
					# would become:
					# IMPRECISE;SVTYPE=INS;MEINFO=L1,-179,180,.,ALU,-16,16,.;TARGETSITEDUPL=duplication,unknown

					# merging of MEINFO and TARGETSITEDUPL in INFO field, eg.
					# if merging 2 of the same type of ALT <INS:ME:L1> then their INFO fields:
					# IMPRECISE;SVTYPE=INS;MEINFO=L1,-16,16,.;TARGETSITEDUPL=duplication
					# IMPRECISE;SVTYPE=INS;MEINFO=L1,-179,180,.;TARGETSITEDUPL=unknown
					# would become the following having the widest region for MEINFO and most information TARGETSITEDUPL:
					# IMPRECISE;SVTYPE=INS;MEINFO=L1,-179,180,.;TARGETSITEDUPL=duplication

					this_alt_fields = this_inline_alt.split(',')
					new_alt_fields_string = previous_ME_inline_alt 			# eg. <INS:ME:ALU> or <INS:ME:ALU>,<INS:ME:L1>
					for i in range( 0, len(this_alt_fields) ):
						this_alt_field = this_alt_fields[i]
						does_this_alt_field_exist_in_list = False
						temp_new_alt_fields = new_alt_fields_string.split(',')
						for j in range( 0, len(temp_new_alt_fields) ):
							temp_new_alt_field = temp_new_alt_fields[j]
							if (temp_new_alt_field == this_alt_field):
								does_this_alt_field_exist_in_list = True
						if (does_this_alt_field_exist_in_list == False):
							new_alt_fields_string = new_alt_fields_string + ',' + this_alt_field

					previous_ME_inline_meinfo = ''
					previous_ME_inline_targetsitedupl = ''
					prev_info_fields = previous_ME_inline_info.split(';')
					for prev_info_field in prev_info_fields:
						bits = prev_info_field.split('=')
						if (bits[0] == 'MEINFO'):
							previous_ME_inline_meinfo = bits[1]
						elif (bits[0] == 'TARGETSITEDUPL'):
							previous_ME_inline_targetsitedupl = bits[1]
					this_inline_meinfo = ''
					this_inline_targetsitedupl = ''
					this_info_fields = this_inline_info.split(';')
					for this_info_field in this_info_fields:
						bits = this_info_field.split('=')
						if (bits[0] == 'MEINFO'):
							this_inline_meinfo = bits[1]
						elif (bits[0] == 'TARGETSITEDUPL'):
							this_inline_targetsitedupl = bits[1]
					# adjust meinfo and targetsitedupl fields to merge best information from the multiple VCF records
					best_meinfo_field = choose_best_meinfo( previous_ME_inline_meinfo, this_inline_meinfo )
					best_targetsitedupl_field = choose_best_targetsitedupl( previous_ME_inline_targetsitedupl, this_inline_targetsitedupl )

					# build the new INFO field with the new versions of MEINFO and TARGETSITEDUPL
					new_info_fields_string = ''
					old_info_fields = previous_ME_inline_info.split(';')
					for old_info_field in old_info_fields:
						bits = old_info_field.split('=')
						if (bits[0] == 'MEINFO'):
							new_info_fields_string = new_info_fields_string + 'MEINFO=' + best_meinfo_field + ';'
						elif (bits[0] == 'TARGETSITEDUPL'):
							new_info_fields_string = new_info_fields_string + 'TARGETSITEDUPL=' + best_targetsitedupl_field + ';'
						else:
							new_info_fields_string = new_info_fields_string + old_info_field + ';'
					if (new_info_fields_string != ''):
						new_len = len(new_info_fields_string) - 1
						new_info_fields_string = new_info_fields_string[0:new_len]
					previous_ME_inline_info = new_info_fields_string

					new_alt_fields = new_alt_fields_string.split(',')
					map_old_alt_to_new_alt = [0] * (len(this_alt_fields) + 1)
					i = 1
					for this_alt_field in this_alt_fields:
						for j in range( 0, len(new_alt_fields) ):
							if (this_alt_field == new_alt_fields[j]):
								new_i = j + 1
								map_old_alt_to_new_alt[i] = new_i
						i = i + 1

					# Now reformat all the formatted samples in this new VCF variant
					# to point to the ALT variants in the new ALT string.
					# Actually, Genotype GT is the only format field to be reformatted.
					# All other format types are just copied across, they don't point to ALT values.
					# Inside the sample's formatted GT field, each value has to be reformatted to point to correct new ALT value.
					# A new string for the formatted samples needs to be produced.
					# It will be the same as the string of formatted samples of the previous inline, except for the following.
					# Where this new VCF variant contains formatted sample for a sample, 
					# that sample's formatted field must come from this new VCF variant record,
					# not from the previous inline.
					# We assume that the GT field is the first field of FORMAT.
					# We will look at the GT fields of the formatted sample field of this new VCF variant record.
					# If it contains any values other than . then we assume that this new VCF variant record contains data for this sample.
					# Otherwise we assume that it does not, and will take the formatted sample field of the previous inline
					# that we are merging this new VCF variant record with.

					# Crash if the first FORMAT field is not GT

					this_inline_format_fields = this_inline_format.split(':')
					if (this_inline_format_fields[0] != 'GT'):
						raise ValueError("\n\nThe first FORMAT field for variant " + str(this_inline_chrom) + ":" + str(this_inline_pos) + " is not GT for Genotype. We are trying to merge multiple ALT variants into one and we need this genotype field to know which VCF record contains each sample's data. Will not continue processing this file.\n")

					# For each sample, figure out whether its data is in the new VCF variant record and thus needs the GT to be reformatted for the merged ALT field,
					# or is its data in the previous inline field so just copy the entire formatted data across from the previous inline field.

					where_is_sample_data_array = []
					previous_inline_sample_data_array = []
					this_inline_multiple_samples_fields_array = this_inline_multiple_samples_fields.split("\t")
					for sample_idx in range( 0, len(this_inline_multiple_samples_fields_array) ):
						previous_ME_inline_multiple_samples_fields_array = previous_ME_inline_multiple_samples_fields.split("\t")
						previous_inline_this_sample = previous_ME_inline_multiple_samples_fields_array[sample_idx]
						this_inline_this_sample = this_inline_multiple_samples_fields_array[sample_idx]
						previous_inline_this_sample_array = previous_inline_this_sample.split(':')
						this_inline_this_sample_array = this_inline_this_sample.split(':')
						previous_inline_this_sample_GT = previous_inline_this_sample_array[0]
						this_inline_this_sample_GT = this_inline_this_sample_array[0]
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT_temp.replace('|','')
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT_temp.replace('/','')
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT_temp.replace('.','')
						where_is_sample_data = 'this_inline'
						if (this_inline_this_sample_GT_temp == ''):
							where_is_sample_data = 'previous_inline'
						where_is_sample_data_array.append( where_is_sample_data )
						previous_inline_sample_data_array.append( previous_inline_this_sample )

					# Now reconstruct the string of formatted sample fields.
					# The formatted data for each sample either comes from the previous_inline, 
					# or from this new VCF inline after being reformatted for the change in its merged ALT alleles

					new_this_inline_multiple_samples_fields_string = ''
					for sample_idx in range( 0, len(this_inline_multiple_samples_fields_array) ):

						formatted_sample_field_for_this_sample = ''
						if (where_is_sample_data_array[sample_idx] == 'previous_inline'):
							formatted_sample_field_for_this_sample = previous_inline_sample_data_array[sample_idx]

						else: # (where_is_sample_data_array[sample_idx] == 'this_inline')

							this_one_sample_multiple_format_fields = this_inline_multiple_samples_fields_array[sample_idx]
							new_one_sample_multiple_format_fields = ''
							this_inline_format_array = this_inline_format.split(':')
							this_inline_one_sample_multiple_format_fields_array = this_one_sample_multiple_format_fields.split(':')
							format_idx = 0
							for this_one_sample_one_format_field in this_inline_one_sample_multiple_format_fields_array:
								new_one_sample_one_format_field = ''
								this_format_type = this_inline_format_array[format_idx]
								if (this_format_type != 'GT'):
									new_one_sample_one_format_field = this_one_sample_one_format_field
								else: # This is the Genotype field, change values according to new ALT positions
									if (this_one_sample_one_format_field == '.'):
										new_one_sample_one_format_field = this_one_sample_one_format_field
									else:
										separator_position = find_first_separator_position( this_one_sample_one_format_field )
										if (separator_position == -1):
											raise ValueError("\n\nSample formatted genotype field for variant " + str(this_inline_chrom) + ":" + str(this_inline_pos) + " does not have a slash / or bar | and we are trying to merge multiple ALT variants into one and accordingly modify this genotype field. Will not continue processing this file.\n")
										sample_GT_fields = re.split( "\||\/", this_one_sample_one_format_field ) # split on | or /
										new_one_sample_one_format_field = ''
										remaining_GT_field_to_look_for_separators = this_one_sample_one_format_field
										for sample_GT_field in sample_GT_fields:
											new_sample_GT_field = '.'
											separator_position, sample_GT_separator = find_first_separator_position( remaining_GT_field_to_look_for_separators )
											remaining_GT_field_to_look_for_separators = remaining_GT_field_to_look_for_separators[ (separator_position+1) : ]
											if (sample_GT_field != '.'):
												new_sample_GT_field = map_old_alt_to_new_alt[int(sample_GT_field)]
												if (new_sample_GT_field == ''):
													raise ValueError("\n\nSample formatted genotype field for variant " + str(this_inline_chrom) + ":" + str(this_inline_pos) + " has a value that can't be matched to one of its ALT variants, and we are trying to modify it to match the new merged ALTs. Will not continue processing this file.\n")
											new_one_sample_one_format_field = new_one_sample_one_format_field + str(new_sample_GT_field) + sample_GT_separator

								# new_one_sample_multiple_format_fields = new_one_sample_multiple_format_fields + ':' + new_one_sample_one_format_field
								if (new_one_sample_one_format_field == ''):
									new_one_sample_one_format_field = '.'
								if (new_one_sample_multiple_format_fields == ''):
									new_one_sample_multiple_format_fields = new_one_sample_one_format_field
								else:
									new_one_sample_multiple_format_fields = new_one_sample_multiple_format_fields + ':' + new_one_sample_one_format_field
								format_idx = format_idx + 1
							formatted_sample_field_for_this_sample = new_one_sample_multiple_format_fields

						# new_this_inline_multiple_samples_fields_string = new_this_inline_multiple_samples_fields_string + "\t" + new_one_sample_multiple_format_fields
						if (new_this_inline_multiple_samples_fields_string == ''):
							new_this_inline_multiple_samples_fields_string = formatted_sample_field_for_this_sample
						else:
							new_this_inline_multiple_samples_fields_string = new_this_inline_multiple_samples_fields_string + "\t" + formatted_sample_field_for_this_sample

					if (new_this_inline_multiple_samples_fields_string == ''):
						new_this_inline_multiple_samples_fields_string = '.'

					previous_ME_inline_alt = new_alt_fields_string
					previous_ME_inline_multiple_samples_fields = new_this_inline_multiple_samples_fields_string

	# We have read and processed all the input VCF records.
	# If there is a previous ME merge group still hanging around then write it out.
	# If there are non-ME inlines still hanging around after them, then write them out.

	if (previous_ME_inline_chrom != ''):

		# Write out the previous ME merge group.
		# It has not been written out yet because we were still looking for more ME variants having similar POS
		# when we hit the end of the VCF file.
		# If necessary get the average position of the list of positions, and it's ref.seq. nucleotide

		previous_ME_inline_avg_pos = previous_ME_inline_pos_array[0]
		if (len(previous_ME_inline_pos_array) > 1):

			# The previous variant is a bunch of variants of similar position.
			# Now that we are writing them out as one merged variant
			# and know all the variants in this merged bunch,
			# get the average POS and corresponding REF

			# Get the average POS

			sum_of_pos = 0
			for this_pos in previous_ME_inline_pos_array:
				sum_of_pos = sum_of_pos + this_pos
			previous_ME_inline_avg_pos = sum_of_pos / len(previous_ME_inline_pos_array)
			previous_ME_inline_avg_pos = int(previous_ME_inline_avg_pos)

			# Get the corresponding REF of that POS

			samtools_faidx_position = str(previous_ME_inline_chrom) + ':' + str(previous_ME_inline_avg_pos) + '-' + str(previous_ME_inline_avg_pos)
			samtools_faidx_command = 'samtools faidx ' + input_reference_genome_fastq_file + ' ' + samtools_faidx_position
			command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
			if (command_status != 0):
				raise ValueError("\n\nWas not able to get REF position from reference genome. Thus will not continue processing any more mobile elements.\n")
			command_output_lines = command_output.split("\n")
			# some positions that mobile element mapped to may not be found in reference genome
			# eg. hs37d5:180728
			if (len(command_output_lines) > 1):
				previous_ME_inline_ref = command_output_lines[1].strip()

		outline = str(previous_ME_inline_chrom) + "\t" + str(previous_ME_inline_avg_pos) + "\t" + str(previous_ME_inline_id) + "\t" + previous_ME_inline_ref + "\t" + previous_ME_inline_alt + "\t" + str(previous_ME_inline_qual) + "\t" + previous_ME_inline_filter + "\t" + previous_ME_inline_info + "\t" + previous_ME_inline_format + "\t" + previous_ME_inline_multiple_samples_fields + "\n"
		sys.stdout.write( outline )

		# Write out any non-ME variants having positions after the previous ME merge group

		if ( len(previous_ME_inline_subsequent_non_ME_inlines_to_write_out) > 0 ):
			for previous_ME_inline_subsequent_non_ME_inline in previous_ME_inline_subsequent_non_ME_inlines_to_write_out:
				sys.stdout.write( previous_ME_inline_subsequent_non_ME_inline )
			previous_ME_inline_subsequent_non_ME_inlines_to_write_out = []

if __name__=='__main__':
    main()


