#!/usr/bin/env python3

import sys
import argparse

def check_locus_span(bedpe_file, locus_chrom, locus_start, locus_end, bases=5000, merge_distance=200000, verbose=False):
    """Check haplotype coverage criteria for a locus."""
    # Read alignments from BEDPE file for specified chromosome
    alignments = []
    contig_names = set()
    try:
        with open(bedpe_file, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 6 and fields[3] == locus_chrom:
                    contig_name = fields[0]
                    contig_names.add(contig_name)
                    alignments.append({
                        'query_name': contig_name,
                        'target_start': int(fields[4]),
                        'target_end': int(fields[5])
                    })
    except FileNotFoundError:
        if verbose:
            print(f"Error: File '{bedpe_file}' not found")
        return False, {}, False
    
    if not alignments:
        if verbose:
            print(f"No alignments found for chromosome {locus_chrom}")
        return False, {}, False
    
    # Sort alignments by start position
    alignments.sort(key=lambda x: x['target_start'])
    
    # CRITERION 1: Check if there's exactly one unique contig
    single_contig = len(contig_names) == 1
    
    # Calculate first and last regions based on fixed number of bases
    first_region_end = locus_start + bases
    last_region_start = locus_end - bases
    
    # CRITERION 2: Check if contigs span first and last regions
    first_covered = False
    last_covered = False
    
    for aln in alignments:
        # Check first region
        if aln['target_start'] <= locus_start and aln['target_end'] >= first_region_end:
            first_covered = True
        
        # Check last region
        if aln['target_start'] <= last_region_start and aln['target_end'] >= locus_end:
            last_covered = True
    
    boundaries_covered = first_covered and last_covered
    
    # CRITERION 3: Check if there's a single contig after merging
    # Merge contigs that are within merge_distance of each other
    merged_contigs = []
    if alignments:
        merged_contigs = [alignments[0].copy()]
        
        for curr in alignments[1:]:
            prev = merged_contigs[-1]
            # If current contig is within merge_distance of previous contig, merge them
            if curr['target_start'] <= prev['target_end'] + merge_distance:
                # Extend the previous contig
                prev['target_end'] = max(prev['target_end'], curr['target_end'])
            else:
                # Add as a new contig
                merged_contigs.append(curr.copy())
    
    single_merged_contig = len(merged_contigs) == 1
    
    # Check if all criteria pass
    all_pass = single_contig and boundaries_covered and single_merged_contig
    
    # Store detailed results
    details = {
        'single_contig': single_contig,
        'boundaries_covered': boundaries_covered,
        'single_merged_contig': single_merged_contig,
        'first_covered': first_covered,
        'last_covered': last_covered
    }
    
    # Print verbose output if requested
    if verbose:
        print(f"\nHaplotype Span Analysis for {locus_chrom}:{locus_start}-{locus_end}")
        print(f"=================================================================")
        print(f"Total alignment regions: {len(alignments)}")
        print(f"Unique contigs: {len(contig_names)}")
        print(f"Regions after merging (d={merge_distance:,}bp): {len(merged_contigs)}")
        
        print(f"\nCRITERIA EVALUATION:")
        print(f"1. Single contig: {'PASS' if single_contig else 'FAIL'}")
        print(f"2. Coverage of first/last {bases:,}bp: {'PASS' if boundaries_covered else 'FAIL'}")
        print(f"   - First {bases:,}bp: {'Covered' if first_covered else 'Not covered'}")
        print(f"   - Last {bases:,}bp: {'Covered' if last_covered else 'Not covered'}")
        print(f"3. Single region after merging: {'PASS' if single_merged_contig else 'FAIL'}")

        print(f"\nOVERALL RESULT: {'PASS' if all_pass else 'FAIL'}")
    
    return all_pass, details, boundaries_covered

def parse_locus(locus_str):
    """Parse a locus string in format chr:start-end."""
    try:
        chrom, positions = locus_str.split(':')
        start, end = map(int, positions.split('-'))
        return (chrom, start, end)
    except:
        raise argparse.ArgumentTypeError("Locus format should be chrom:start-end (e.g., chr6:31972057-32055418)")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Check if alignments in a BEDPE file meet coverage criteria for a locus.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -b alignments.bedpe -l chr6:31972057-32055418
  %(prog)s -v -b alignments.bedpe -l chr6:31972057-32055418 -bp 1000 -d 100000
"""
    )
    
    # Required arguments with flags
    parser.add_argument('-b', '--bed', dest='bedpe_file', required=True,
                      help='Path to the BEDPE file with alignments')
    parser.add_argument('-l', '--locus', required=True, type=parse_locus,
                      help='Reference locus in format chr:start-end')
    
    # Optional arguments
    parser.add_argument('-v', '--verbose', action='store_true', 
                      help='Show detailed report')
    parser.add_argument('-bp', '--bases', type=int, default=1000,
                      help='Number of base pairs to check at start and end of locus (must be > 0 and < 50%% of locus length, default: 1000)')
    parser.add_argument('-d', '--distance', type=int, default=200000,
                      help='Merge distance in bp (default: 200000)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Extract locus components
    chrom, start, end = args.locus
    
    # Validate the bases parameter after we know the locus length
    locus_length = end - start
    max_bases = locus_length // 2  # 50% of locus length
    
    if args.bases <= 0 or args.bases >= max_bases:
        parser.error(f"Bases must be > 0 and < 50% of locus length ({max_bases})")

    
    # Run the check
    all_pass, details, boundaries_pass = check_locus_span(
        args.bedpe_file, chrom, start, end, 
        args.bases, args.distance, args.verbose)
    
    # Output true/false to stdout for scripting
    if not args.verbose:
        print("true" if all_pass else "false")
    
    # Use exit code for programmatic checking
    sys.exit(0 if all_pass else 1)
