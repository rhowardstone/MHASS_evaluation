#!/usr/bin/env python3

import argparse
import os
import sys
import math
from pathlib import Path


def read_fasta_generator(fasta_file):
    """
    Generator function to read a FASTA file entry by entry.
    Yields tuples of (sequence_id, sequence).
    """
    with open(fasta_file, 'r') as f:
        current_id = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Yield the previous sequence before starting a new one
                if current_id:
                    yield (current_id, ''.join(current_seq))
                
                current_id = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget to yield the last sequence
        if current_id:
            yield (current_id, ''.join(current_seq))


def count_sequences(fasta_file):
    """Count the number of sequences in a FASTA file."""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def split_fasta(fasta_file, output_dir, prefix, num_chunks=None, seqs_per_chunk=None):
    """
    Split a FASTA file into multiple chunks.
    
    Parameters:
    - fasta_file: Path to the input FASTA file
    - output_dir: Directory where output chunks will be saved
    - prefix: Prefix for output filenames
    - num_chunks: Number of chunks to create (mutually exclusive with seqs_per_chunk)
    - seqs_per_chunk: Number of sequences per chunk (mutually exclusive with num_chunks)
    
    Returns:
    - List of output filenames
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Count sequences if needed
    if num_chunks and not seqs_per_chunk:
        total_seqs = count_sequences(fasta_file)
        seqs_per_chunk = math.ceil(total_seqs / num_chunks)
        print(f"Total sequences: {total_seqs}, sequences per chunk: {seqs_per_chunk}")
    
    # Initialize variables
    current_chunk = 1
    current_seq_count = 0
    output_files = []
    current_outfile = None
    
    # Process the FASTA file
    for seq_id, sequence in read_fasta_generator(fasta_file):
        # Open a new file if needed
        if current_seq_count == 0:
            chunk_filename = os.path.join(output_dir, f"{prefix}_{current_chunk:04d}.fa")
            output_files.append(chunk_filename)
            current_outfile = open(chunk_filename, 'w')
        
        # Write the sequence to the current chunk
        current_outfile.write(f"{seq_id}\n{sequence}\n")
        current_seq_count += 1
        
        # Close current file and reset counter if chunk is complete
        if current_seq_count >= seqs_per_chunk:
            current_outfile.close()
            current_chunk += 1
            current_seq_count = 0
    
    # Close the last file if it's still open
    if current_outfile and not current_outfile.closed:
        current_outfile.close()
    
    return output_files


def main():
    parser = argparse.ArgumentParser(description='Split a FASTA file into multiple chunks')
    parser.add_argument('fasta_file', help='Input FASTA file to split')
    parser.add_argument('-o', '--output-dir', default='./chunks', help='Directory for output chunks (default: ./chunks)')
    parser.add_argument('-p', '--prefix', default='chunk', help='Prefix for output filenames (default: chunk)')
    
    # Mutually exclusive options for chunk size
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-n', '--num-chunks', type=int, help='Number of chunks to create')
    group.add_argument('-s', '--seqs-per-chunk', type=int, help='Number of sequences per chunk')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.fasta_file):
        print(f"Error: Input file '{args.fasta_file}' does not exist", file=sys.stderr)
        sys.exit(1)
        
    print(f"Splitting FASTA file: {args.fasta_file}")
    
    # Split the file
    output_files = split_fasta(
        args.fasta_file, 
        args.output_dir, 
        args.prefix, 
        num_chunks=args.num_chunks, 
        seqs_per_chunk=args.seqs_per_chunk
    )
    
    # Print summary
    print(f"\nSplit complete!")
    print(f"Created {len(output_files)} chunks in directory: {args.output_dir}")
    #print(f"Files:")
    #for i, f in enumerate(output_files, 1):
    #    print(f"  {i}. {os.path.basename(f)}")


if __name__ == "__main__":
    main()
