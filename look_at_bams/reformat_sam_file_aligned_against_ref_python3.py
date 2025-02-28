#!/usr/bin/env python3

import sys

######################################################
def is_integer(s):
    """Check if a given string is an integer."""
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def read_fasta_file(input_fasta_ref_file_name):
    """Reads a FASTA reference file and extracts metadata and sequence."""
    with open(input_fasta_ref_file_name, 'r') as input_fasta_ref:
        input_fasta_ref_lines = input_fasta_ref.readlines()

    fasta_hdr = input_fasta_ref_lines[0][1:].replace('-', ':').strip()
    fasta_hdr_fields = fasta_hdr.split(':')

    ref_name = 'human'
    ref_chrom = fasta_hdr_fields[0]
    ref_start = int(fasta_hdr_fields[1])
    ref_end = int(fasta_hdr_fields[2])

    ref_sequence = ''.join(line.strip() for line in input_fasta_ref_lines[1:])

    return ref_name, ref_start, ref_end, ref_sequence

######################################################
def process_sam_file(input_sam_file_name):
    """Reads a SAM file and extracts relevant information."""
    with open(input_sam_file_name, 'r') as input_sam:
        return [line.strip() for line in input_sam.readlines() if not line.startswith('@')]

######################################################
def generate_reference_display(ref_start, ref_end, ref_name, ref_sequence):
    """Generates formatted reference sequence for alignment."""
    reference_array = []

    # Initialize the reference array
    reference_array = []

    # First Line: Reference Positions
    ref_start_display = f"{ref_start}:"
    leftcol = f"{ref_start_display:>130}"  # Right-align ref_start
    grid = ""

    for i in range(ref_start, ref_end + 1):
        grid_char = " "

        # Remove markers at 7, 8, 9 for proper spacing
        if i % 10 in {7, 8, 9}:
            grid_char = ""

        # Insert numbers at multiples of 10, keeping only the last 4 digits
        if i % 10 == 0:
            grid_char = str(i)[-4:]

            # Adjust for first few numbers at the start of the reference
            if i <= (ref_start + 2):
                chr_len = len(grid_char)
                if i == ref_start:
                    grid_char = grid_char[-1:]
                elif i == (ref_start + 1):
                    grid_char = grid_char[-2:]
                elif i == (ref_start + 2):
                    grid_char = grid_char[-3:]

        grid += grid_char

    outline = f"{leftcol}{grid}:{ref_end} \n"
    reference_array.append(outline)

    # Second Line: Grid Markers
    leftcol = f"{ref_name:<130}"  # Left-align reference name
    grid = ""

    for i in range(ref_start, ref_end + 1):
        grid_char = "-"

        # Place '+' at multiples of 5
        if i % 5 == 0:
            grid_char = "+"

        # Place numbers at multiples of 10, keeping only the second-last digit
        if i % 10 == 0:
            grid_char = str(i)[-2:-1]

        grid += grid_char

    outline = f"{leftcol}{grid}\n"
    reference_array.append(outline)

    # Third Line: Reference Sequence
    leftcol = f"{ref_name:<130}"
    outline = f"{leftcol}{ref_sequence}\n"
    reference_array.append(outline)

    return reference_array

######################################################
def process_sam_reads(input_sam_lines, ref_start, reference_array):
    """Processes SAM file reads and aligns them under the reference sequence."""
    outline_array = list(reference_array)

    # Output the SAM read alignments
    seq_line_count = 0

    for inline in input_sam_lines:
        inline = inline.strip()

        # Ignore header lines
        if not inline.startswith('@'):
            infields = inline.split("\t")

            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = infields[:11]

            # Format each field to match original alignment
            qname_pr = f"{qname:<54}"  # Left-align in a 54-character space
            rname_pr = f"{rname:<4}"   # Left-align in a 4-character space
            pos_pr = f"{pos:<12}"      # Left-align in a 12-character space
            cigar_pr = f"{cigar:<12}"  # Left-align in a 12-character space
            rnext_pr = f"{rnext:<4}"   # Left-align in a 4-character space
            pnext_pr = f"{pnext:<12}"  # Left-align in a 12-character space

            # Ensure CIGAR string does not exceed 12 characters
            if len(cigar_pr) > 12:
                cigar_pr = cigar_pr[:12]

            # Construct left column
            leftcol = f"{qname_pr}{rname_pr}{pos_pr}{cigar_pr}{rnext_pr}{pnext_pr}"
            leftcol = f"{leftcol:<130}"  # Ensure proper left alignment

            # Determine soft-clipped bases at start of read
            first_soft_clip = 0
            cigar_split = cigar.split('S')

            if is_integer(cigar_split[0]):
                first_soft_clip = int(cigar_split[0])

            # Calculate padding required for correct alignment
            padding_needed = int(pos) - ref_start - first_soft_clip
            padding_needed = max(0, padding_needed)  # Ensure padding is non-negative

            padding = ' ' * padding_needed
            outline = f"{leftcol}{padding}{seq}\n"
            outline_array.append(outline)

            # Output the reference sequence every 50 lines of SAM fragments
            seq_line_count += 1
            if seq_line_count % 50 == 0:
                for refline in reference_array:
                    outline_array.append(refline)
                seq_line_count = 0

    return outline_array

######################################################
def main():
    """Main function to execute the script logic."""
    if len(sys.argv) != 3:
        print("Usage: python3 reformat_sam_file_aligned_against_ref.py <input_sam> <input_fasta_ref>")
        sys.exit(1)

    input_sam_file_name = sys.argv[1].strip()
    input_fasta_ref_file_name = sys.argv[2].strip()

    # Read input FASTA reference
    ref_name, ref_start, ref_end, ref_sequence = read_fasta_file(input_fasta_ref_file_name)

    # Read input SAM file
    input_sam_lines = process_sam_file(input_sam_file_name)

    # Generate reference sequence alignment
    reference_array = generate_reference_display(ref_start, ref_end, ref_name, ref_sequence)

    # Process SAM file reads
    outline_array = process_sam_reads(input_sam_lines, ref_start, reference_array)

    # Write to output file
    outfile_name = f"{input_sam_file_name}.aligned.txt"
    with open(outfile_name, 'w') as outfile:
        outfile.writelines(outline_array)
        outfile.write("\n\n")

    print(f"Alignment saved to {outfile_name}")

######################################################
if __name__ == '__main__':
    main()

