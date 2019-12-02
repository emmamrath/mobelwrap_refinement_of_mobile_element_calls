#!/usr/bin/python
# python reformat_sam_file_aligned_against_ref_first_get_extracts.py what_is_the_extract where_is_bam where_is_ref_fasta
# python reformat_sam_file_aligned_against_ref_first_get_extracts.py 1:10225 /nvme/mobile_elements_2016dec/debug_MGRB_picard_header/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.bam /nvme/mobile_elements_2016dec/1000_genomes_hs37d5_fastq/human_g1k_v37.fasta
# python reformat_sam_file_aligned_against_ref_first_get_extracts.py 1:564336:564704 /nvme/mobile_elements_2016dec/debug_MGRB_picard_header/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.bam /nvme/mobile_elements_2016dec/1000_genomes_hs37d5_fastq/human_g1k_v37.fasta
# python reformat_sam_file_aligned_against_ref_first_get_extracts.py 1:564336-564704 /nvme/mobile_elements_2016dec/debug_MGRB_picard_header/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.bam /nvme/mobile_elements_2016dec/1000_genomes_hs37d5_fastq/human_g1k_v37.fasta

# /g/data3/abc/results/phase2/hs37d5x/GATK3_fastq2gvcf-hs37d5x-1.0/bams/

# This program generates the following commands and runs these commands to produce files:
# 	samtools view -h /nvme/mobile_elements_2016dec/debug_MGRB_picard_header/MGRBp1_sample1000_sg1_humansingle1_999.picard_sorted.dupmarked.bam 1:900-1100 > MGRBp1_sample1000_extract_1_900_1100.sam
# 	samtools faidx human_g1k_v37.fasta 1:700-1300 > human_g1k_v37_1_700_1300.fasta
# 	python reformat_sam_file_aligned_against_ref.py MGRBp1_sample1000_extract_1_900_1100.sam human_g1k_v37_1_700_1300.fasta

# The output of the last program is the bam fragments visually aligned to the reference fasta.


import sys
import os

def main():

	### Read input parameters

	what_is_the_extract = sys.argv[1]
	where_is_bam = sys.argv[2]
	where_is_ref_fasta = sys.argv[3]
	what_is_the_extract = what_is_the_extract.strip()
	where_is_bam = where_is_bam.strip()
	where_is_ref_fasta = where_is_ref_fasta.strip()

	### Read input parameters

	what_is_the_extract = what_is_the_extract.replace('_',':')
	what_is_the_extract = what_is_the_extract.replace('-',':')
	infields = what_is_the_extract.split(':')
	chrom = infields[0]
	# chrom = chrom.replace( 'chr', '' )
	num_in = int(infields[1])
	num_in_2 = num_in
	if (len(infields) >= 3):
		num_in_2 = int(infields[2])
	num1 = num_in - 100
	num2 = num_in_2 + 100
	num3 = num_in - 300
	num4 = num_in_2 + 300
	if (num1 < 1):
		num1 = 1
	if (num3 < 1):
		num3 = 1
	num1 = str(num1)
	num2 = str(num2)
	num3 = str(num3)
	num4 = str(num4)

	# generate the commands to run

	command1 = 'samtools view ' + where_is_bam + ' ' + chrom + ':' + num1 + '-' + num2 + ' > temp_bam_extract_' + chrom + '_' + num1 + '_' + num2 + '.sam'
	command2 = 'samtools faidx ' + where_is_ref_fasta + ' ' + chrom + ':' + num3 + '-' + num4 + ' > temp_fasta_extract_' + chrom + '_' + num3 + '_' + num4 + '.fasta'
	command3 = 'python reformat_sam_file_aligned_against_ref.py temp_bam_extract_' + chrom + '_' + num1 + '_' + num2 + '.sam temp_fasta_extract_' + chrom + '_' + num3 + '_' + num4 + '.fasta'

	print 'Here are the 3 commands that will now be run:'
	print command1
	print command2
	print command3

	# run the commands

	os.system( command1 )
	os.system( command2 )
	os.system( command3 )


if __name__=='__main__':
    main()

