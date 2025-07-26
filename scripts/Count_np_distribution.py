#!/usr/bin/env python3
"""
Extract number of passes (NP) and read length from PacBio HiFi FASTQ files in parallel.

PacBio HiFi read headers contain the number of passes in the format:
@movie/zmw/ccs np=X

This script extracts NP and read length for all reads across multiple FASTQ files.
"""

import os
import gzip
import argparse
from multiprocessing import Pool, cpu_count
from pathlib import Path
import re
import sys

def process_fastq_file(args_tuple):
    """
    Process a single FASTQ file and extract (np, read_length) pairs.
    Filters: NP < max_np (if specified) and read length between user-defined minlen and maxlen.
    
    Returns:
        list of tuples: [(np, read_length), ...]
    """
    filepath, minlen, maxlen, max_np = args_tuple
    results = []
    total_reads = 0
    filtered_reads = 0
    
    # Pattern to extract np from header: np=(\d+)
    np_pattern = re.compile(r'np=(\d+)')
    
    # Determine if file is gzipped
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    try:
        with open_func(filepath, mode) as f:
            line_count = 0
            header = None
            
            for line in f:
                line_count += 1
                line = line.strip()
                
                # FASTQ format: header, sequence, +, quality
                remainder = (line_count - 1) % 4
                
                if remainder == 0:  # Header line
                    header = line
                elif remainder == 1:  # Sequence line
                    if header and header.startswith('@'):
                        # Extract NP from header
                        match = np_pattern.search(header)
                        if match:
                            np = int(match.group(1))
                            read_length = len(line)
                            total_reads += 1
                            
                            # Apply filters: NP < max_np (if specified) and length in [minlen, maxlen]
                            np_passes = max_np is None or np < max_np
                            length_passes = minlen <= read_length <= maxlen
                            
                            if np_passes and length_passes:
                                results.append((np, read_length))
                                filtered_reads += 1
                        else:
                            print(f"Warning: No 'np=' found in header: {header}", file=sys.stderr)
        
        print(f"Processed {filepath}: {total_reads} total reads, {filtered_reads} passed filters")
        return results
    
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)
        return []

def find_fastq_files(directory):
    """Find all FASTQ files in the given directory."""
    fastq_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    fastq_files = []
    
    for ext in fastq_extensions:
        fastq_files.extend(Path(directory).glob(f'*{ext}'))
    
    return [str(f) for f in fastq_files]

def main():
    parser = argparse.ArgumentParser(
        description='Extract NP and read length from PacBio HiFi FASTQ files'
    )
    parser.add_argument('input_dir', help='Directory containing FASTQ files')
    parser.add_argument('output', help='Output TSV file with np, read_length columns')
    parser.add_argument('-j', '--jobs', type=int, default=None,
                       help='Number of parallel jobs (default: number of CPUs)')
    parser.add_argument('--pattern', default='*',
                       help='File pattern to match (default: * for all FASTQ files)')
    parser.add_argument('--np-summary', help='Output file for NP count summary (optional)')
    parser.add_argument('--minlen', type=int, default=2000,
                       help='Minimum read length to include (default: 2000)')
    parser.add_argument('--maxlen', type=int, default=3000,
                       help='Maximum read length to include (default: 3000)')
    parser.add_argument('--max-np', type=int, default=None,
                       help='Maximum number of passes to include (default: no filtering)')
    
    args = parser.parse_args()
    
    # Find FASTQ files
    print(f"Searching for FASTQ files in {args.input_dir}...")
    if args.pattern != '*':
        fastq_files = list(Path(args.input_dir).glob(args.pattern))
        fastq_files = [str(f) for f in fastq_files if f.suffix in ['.fastq', '.fq', '.gz']]
    else:
        fastq_files = find_fastq_files(args.input_dir)
    
    if not fastq_files:
        print("Error: No FASTQ files found!", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(fastq_files)} FASTQ files")
    
    # Display filter information
    filter_desc = f"read length in [{args.minlen}, {args.maxlen}] bp"
    if args.max_np is not None:
        filter_desc = f"NP < {args.max_np} and " + filter_desc
    print(f"Applying filters: {filter_desc}")
    
    # Process files in parallel
    n_jobs = args.jobs or cpu_count()
    print(f"Processing with {n_jobs} parallel jobs...")
    
    input_args = [(f, args.minlen, args.maxlen, args.max_np) for f in fastq_files]
    with Pool(n_jobs) as pool:
        all_results = pool.map(process_fastq_file, input_args)
    
    # Combine results
    print("\nCombining results...")
    combined_results = []
    for results in all_results:
        combined_results.extend(results)
    
    print(f"Total reads passing filters: {len(combined_results)}")
    
    # Write detailed output
    print(f"Writing detailed results to {args.output}...")
    with open(args.output, 'w') as f:
        f.write("np\tread_length\n")
        for np, read_length in combined_results:
            f.write(f"{np}\t{read_length}\n")
    
    # Create NP summary if requested
    if args.np_summary:
        print(f"Creating NP count summary in {args.np_summary}...")
        np_counts = {}
        for np, _ in combined_results:
            np_counts[np] = np_counts.get(np, 0) + 1
        
        with open(args.np_summary, 'w') as f:
            for np in sorted(np_counts.keys()):
                f.write(f"{np}\t{np_counts[np]}\n")
        
        print(f"NP summary written to {args.np_summary}")
    
    # Print summary statistics
    if combined_results:
        nps = [r[0] for r in combined_results]
        lengths = [r[1] for r in combined_results]
        
        print("\nSummary statistics (filtered data):")
        print(f"NP range: {min(nps)} - {max(nps)}")
        print(f"Read length range: {min(lengths)} - {max(lengths)}")
        print(f"Mean NP: {sum(nps)/len(nps):.2f}")
        print(f"Mean read length: {sum(lengths)/len(lengths):.0f}")
        
        if len(set(lengths)) > 1:
            from statistics import correlation
            try:
                corr = correlation(nps, lengths)
                print(f"Correlation (NP, length): {corr:.3f}")
            except:
                print("Could not compute correlation")
    
    print("\nDone!")

if __name__ == "__main__":
    main()
