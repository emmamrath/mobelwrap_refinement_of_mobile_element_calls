#!/usr/bin/env python3

import sys
import os

def main():
    """ Main function to extract BAM and FASTA regions, then reformat the SAM file. """

    # Read input parameters
    if len(sys.argv) != 4:
        print("Usage: python3 reformat_sam_file_aligned_against_ref_first_get_extracts.py <extract_region> <bam_file> <reference_fasta>")
        sys.exit(1)

    extract_region = sys.argv[1].strip()
    bam_file = sys.argv[2].strip()
    ref_fasta = sys.argv[3].strip()

    # Normalize extract region formatting
    extract_region = extract_region.replace('_', ':').replace('-', ':')
    infields = extract_region.split(':')

    if len(infields) < 2:
        print("Error: Invalid extract region format. Expected format: chrom:start or chrom:start-end")
        sys.exit(1)

    chrom = infields[0]
    num_in = int(infields[1])
    num_in_2 = num_in if len(infields) < 3 else int(infields[2])

    # Define extraction ranges
    num1 = max(1, num_in - 100)
    num2 = num_in_2 + 100
    num3 = max(1, num_in - 300)
    num4 = num_in_2 + 300

    # Convert to string for use in commands
    num1, num2, num3, num4 = map(str, (num1, num2, num3, num4))

    # Generate the shell commands
    sam_extract = f"temp_bam_extract_{chrom}_{num1}_{num2}.sam"
    fasta_extract = f"temp_fasta_extract_{chrom}_{num3}_{num4}.fasta"

    command1 = f"samtools view {bam_file} {chrom}:{num1}-{num2} > {sam_extract}"
    command2 = f"samtools faidx {ref_fasta} {chrom}:{num3}-{num4} > {fasta_extract}"
    command3 = f"python3 reformat_sam_file_aligned_against_ref_python3.py {sam_extract} {fasta_extract}"

    # Display commands
    print("Here are the 3 commands that will now be run:")
    print(command1)
    print(command2)
    print(command3)

    # Run the commands
    os.system(command1)
    os.system(command2)
    os.system(command3)

if __name__ == '__main__':
    main()

