#!/usr/bin/env python3
"""
bam_np_length.py
Extract number of passes (NP) and read length from PacBio HiFi CCS BAM/CRAM files
in parallel, mirroring the behaviour of the existing FASTQ-based script.
Header tag expected:  np:i:<int>
Outputs
  1. A detailed TSV  (np<TAB>read_length)
  2. Optional NP histogram  (same layout as original script)
Usage example
  python bam_np_length.py  ccs_dir/  results.tsv \
      --np-summary np_counts.tsv  --minlen 2000 --maxlen 3000 --max-np 60 -j 8
"""
import argparse
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
import pysam   # pip install pysam

def process_bam_file(args_tuple):
    """
    Scan one BAM/CRAM, returning (np, read_length) tuples that pass filters.
    """
    filepath, minlen, maxlen, max_np = args_tuple
    results = []
    total, kept = 0, 0
    
    try:
        with pysam.AlignmentFile(filepath, "rb", check_sq=False) as bam:
            for read in bam:
                total += 1
                if not read.has_tag("np"):
                    continue
                
                np_val = read.get_tag("np")
                length = read.query_length
                
                # Apply filters: NP < max_np (if specified) and length in [minlen, maxlen]
                np_passes = max_np is None or np_val < max_np
                length_passes = minlen <= length <= maxlen
                
                if np_passes and length_passes:
                    results.append((np_val, length))
                    kept += 1
                    
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)
        return []
    
    print(f"Processed {filepath}: {total} reads, {kept} passed filters")
    return results

def find_bam_files(directory, pattern):
    """Return a list of BAM/CRAM files in *directory* matching *pattern*."""
    exts = (".bam", ".cram")
    paths = Path(directory).glob(pattern)
    return [str(p) for p in paths if p.suffix in exts]

def main():
    p = argparse.ArgumentParser(
        description="Extract NP and read length from PacBio HiFi BAM/CRAM files"
    )
    p.add_argument("input_dir", help="Directory containing BAM/CRAM files")
    p.add_argument("output", help="Output TSV file (np, read_length)")
    p.add_argument("-j", "--jobs", type=int, default=None,
                   help="Parallel jobs (default: #CPUs)")
    p.add_argument("--pattern", default="*",
                   help="Filename glob pattern (default: *)")
    p.add_argument("--np-summary", help="Write NP histogram to this file")
    p.add_argument("--minlen", type=int, default=2000,
                   help="Minimum read length (default 2000)")
    p.add_argument("--maxlen", type=int, default=3000,
                   help="Maximum read length (default 3000)")
    p.add_argument("--max-np", type=int, default=None,
                   help="Maximum number of passes to include (default: no filtering)")
    
    args = p.parse_args()
    
    bam_files = find_bam_files(args.input_dir, args.pattern)
    if not bam_files:
        sys.exit("No BAM/CRAM files found")
    
    print(f"Found {len(bam_files)} BAM/CRAM files")
    
    # Display filter information
    filter_desc = f"length {args.minlen}–{args.maxlen}"
    if args.max_np is not None:
        filter_desc = f"NP < {args.max_np}, " + filter_desc
    print(f"Filters: {filter_desc}")
    
    n_jobs = args.jobs or cpu_count()
    with Pool(n_jobs) as pool:
        combined = []
        for res in pool.map(process_bam_file,
                            [(f, args.minlen, args.maxlen, args.max_np) for f in bam_files]):
            combined.extend(res)
    
    print(f"Total reads passing filters: {len(combined)}")
    
    # detailed TSV
    with open(args.output, "w") as out:
        out.write("np\tread_length\n")
        for np_val, length in combined:
            out.write(f"{np_val}\t{length}\n")
    print(f"Wrote detailed results to {args.output}")
    
    # optional NP histogram
    if args.np_summary:
        counts = {}
        for np_val, _ in combined:
            counts[np_val] = counts.get(np_val, 0) + 1
        with open(args.np_summary, "w") as out:
            for np_val in sorted(counts):
                out.write(f"{np_val}\t{counts[np_val]}\n")
        print(f"Wrote NP summary to {args.np_summary}")

if __name__ == "__main__":
    main()
