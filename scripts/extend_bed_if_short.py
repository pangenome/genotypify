#!/usr/bin/env python3

import argparse
import sys

def load_fasta_index(fai_file):
    """Load chromosome lengths from a FASTA index file."""
    chr_lengths = {}
    try:
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chr_name = parts[0]
                    chr_len = int(parts[1])
                    chr_lengths[chr_name] = chr_len
    except Exception as e:
        sys.stderr.write(f"Error reading FASTA index file: {e}\n")
        sys.exit(1)
    return chr_lengths

def extend_bed_regions(bed_stream, chr_lengths, min_len, output_stream):
    """Extend BED regions to satisfy minimum length requirements."""
    for line in bed_stream:
        # Skip comment or header lines
        if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
            output_stream.write(line)
            continue

        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            # Not a valid BED line, write as is
            output_stream.write(line)
            continue

        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])

        # If the chromosome isn't in our index, we can't extend it safely
        if chrom not in chr_lengths:
            output_stream.write(line)
            continue

        chr_len = chr_lengths[chrom]
        region_len = end - start

        # If the region is already long enough, no need to extend
        if region_len >= min_len:
            output_stream.write(line)
            continue

        # Calculate how much to extend
        extend_total = min_len - region_len
        extend_each_side = extend_total // 2
        extra_base = extend_total % 2  # Handle odd-length extensions

        # Calculate potential new boundaries
        new_start = max(0, start - extend_each_side)
        new_end = min(chr_len, end + extend_each_side + extra_base)  # Add extra base to end

        # Check if we hit chromosome boundaries and need to extend more on the other side
        start_extension = start - new_start
        end_extension = new_end - end

        # If we couldn't extend fully on the start side, add more to the end
        if start_extension < extend_each_side and new_end < chr_len:
            extra_needed = extend_each_side - start_extension
            new_end = min(chr_len, new_end + extra_needed)

        # If we couldn't extend fully on the end side, add more to the start
        if end_extension < extend_each_side + extra_base and new_start > 0:
            extra_needed = (extend_each_side + extra_base) - end_extension
            new_start = max(0, new_start - extra_needed)

        # Final check to ensure we've reached minimum length
        final_length = new_end - new_start
        if final_length < min_len and new_end < chr_len:
            new_end = min(chr_len, new_start + min_len)

        # Update the BED line with new coordinates
        parts[1] = str(new_start)
        parts[2] = str(new_end)
        output_stream.write('\t'.join(parts) + '\n')

def main():
    parser = argparse.ArgumentParser(
        description='Extend BED regions (read from stdin) to satisfy minimum length requirements.'
    )
    parser.add_argument(
        'fai_file',
        help='FASTA index file (.fai) for chromosome lengths'
    )
    parser.add_argument(
        'min_len',
        type=int,
        help='Minimum length for regions'
    )
    args = parser.parse_args()

    # Load chromosome lengths
    chr_lengths = load_fasta_index(args.fai_file)

    # Read BED from stdin, write extended BED to stdout
    extend_bed_regions(sys.stdin, chr_lengths, args.min_len, sys.stdout)

if __name__ == "__main__":
    main()
