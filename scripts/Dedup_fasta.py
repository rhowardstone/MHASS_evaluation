#!/usr/bin/env python3
import argparse
import collections
import re
import os
from typing import Dict, List, Tuple, Set

def parse_args():
    parser = argparse.ArgumentParser(description='De-duplicate sequences in a FASTA file with genome-aware grouping.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output de-duplicated FASTA file')
    parser.add_argument('-m', '--mapping', required=True, help='Output mapping file (dedup_id -> original_ids)')
    parser.add_argument('-g', '--genome_labels', help='Optional genome labels file for genome-aware grouping')
    parser.add_argument('-s', '--separator', default='/', help='Separator for joined headers (default: "/")')
    return parser.parse_args()

def extract_source_filename(header: str) -> str:
    """Extract the source filename from a sequence header."""
    # Split by '.source=' and take the second part
    if '.source=' in header:
        parts = header.split('.source=')
        if len(parts) > 1:
            # Now split by '.coordinates=' and take the first part
            source_part = parts[1].split('.coordinates=')[0]
            return source_part
    return None

def extract_coordinates(header: str) -> str:
    """Extract coordinates from header for sorting."""
    match = re.search(r'coordinates=(\d+)-(\d+)', header)
    if match:
        return f"{int(match.group(1)):09d}-{int(match.group(2)):09d}"
    return header

def create_clean_id(headers: List[str], genome: str, dedup_index: int) -> str:
    """
    Create a clean ID using the format: {genome}-{indices}
    Where indices are joined with - for multi-copy ASVs.
    For example: 
    - Single copy: Escherichia_coli-5
    - Multi-copy with consecutive indices: Escherichia_coli-1-2-3
    """
    # If there are multiple headers (indicating identical sequences from different locations),
    # the ID should include all indices separated by dashes
    if len(headers) > 1:
        # Use a range of consecutive indices
        start_idx = dedup_index
        end_idx = start_idx + len(headers) - 1
        indices = '-'.join(str(i) for i in range(start_idx, end_idx + 1))
        return f"{genome}-{indices}"
    else:
        # Single index for single copy
        return f"{genome}-{dedup_index}"

def read_genome_labels(labels_file: str) -> Dict[str, str]:
    """Read genome labels from file, mapping source filenames to genome IDs."""
    labels = {}
    if not labels_file:
        return labels
    
    print(f"Reading genome labels from {labels_file}")
    with open(labels_file, 'r') as f:
        # Skip header if present
        first_line = f.readline().strip()
        if not first_line.startswith('asvid') and '\t' in first_line:
            # Process first line if it's not a header
            parts = first_line.split('\t')
            if len(parts) >= 2:
                labels[parts[0]] = parts[1]
        
        # Process remaining lines
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                source_file = parts[0]  # This is the filename
                genome_id = parts[1]    # This is the genome ID we want to use
                labels[source_file] = genome_id
    
    print(f"Loaded {len(labels)} genome labels")
    return labels

def dedup_fasta_genome_aware(input_file: str, output_file: str, mapping_file: str, 
                             genome_labels: Dict[str, str] = None, separator: str = '/') -> None:
    """De-duplicate FASTA sequences with genome-aware grouping."""
    
    # Dictionary to store sequences and their headers, grouped by genome
    genome_seq_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    
    # Parse input FASTA
    current_header = None
    current_seq = []
    
    print(f"Reading sequences from {input_file}")
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Store previous sequence if it exists
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    
                    # Extract source filename from header
                    source_file = extract_source_filename(current_header)
                    
                    # Look up genome ID from labels
                    if genome_labels and source_file in genome_labels:
                        genome = genome_labels[source_file]
                    else:
                        # Fallback: use source filename if no labels provided
                        genome = source_file if source_file else "Unknown"
                        print(f"Warning: No genome label found for source: {source_file}")
                    
                    genome_seq_dict[genome][seq].append(current_header)
                
                # Start new sequence
                current_header = line[1:]  # Remove the '>' character
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget to store the last sequence
        if current_header and current_seq:
            seq = ''.join(current_seq)
            
            # Extract source filename from header
            source_file = extract_source_filename(current_header)
            
            # Look up genome ID from labels
            if genome_labels and source_file in genome_labels:
                genome = genome_labels[source_file]
            else:
                genome = source_file if source_file else "Unknown"
                print(f"Warning: No genome label found for source: {source_file}")
            
            genome_seq_dict[genome][seq].append(current_header)
    
    # Create mapping data
    mapping_data = {}  # clean_id -> (genome, original_headers, sequence)
    
    # Process sequences by genome
    print("De-duplicating sequences by genome...")
    for genome in sorted(genome_seq_dict.keys()):
        seq_dict = genome_seq_dict[genome]
        print(f"  {genome}: {len(seq_dict)} unique sequences")
        
        # Sort sequences by coordinates for consistent ordering
        sorted_seqs = sorted(seq_dict.items(), 
                            key=lambda x: extract_coordinates(x[1][0]))
        
        # Track the current index for each genome
        current_idx = 1
        
        for seq, headers in sorted_seqs:
            clean_id = create_clean_id(headers, genome, current_idx)
            mapping_data[clean_id] = (genome, headers, seq)
            
            # Increment the index based on the number of copies (headers)
            if len(headers) > 1:
                current_idx += len(headers)
            else:
                current_idx += 1
    
    # Write de-duplicated sequences to output file
    print(f"Writing {len(mapping_data)} de-duplicated sequences to {output_file}")
    with open(output_file, 'w') as f:
        for clean_id in sorted(mapping_data.keys()):
            genome, headers, seq = mapping_data[clean_id]
            f.write(f'>{clean_id}\n')
            
            # Write sequence with line breaks every 60 characters
            for i in range(0, len(seq), 60):
                f.write(f'{seq[i:i+60]}\n')
    
    # Write mapping file
    print(f"Writing ID mapping to {mapping_file}")
    with open(mapping_file, 'w') as f:
        f.write("dedup_id\tgenome\tnum_original\toriginal_headers\n")
        for clean_id in sorted(mapping_data.keys()):
            genome, headers, seq = mapping_data[clean_id]
            # Join headers with separator
            headers_str = separator.join(headers)
            f.write(f"{clean_id}\t{genome}\t{len(headers)}\t{headers_str}\n")
    
    # Print summary statistics
    print("\nSummary:")
    total_original = sum(len(headers) for _, (_, headers, _) in mapping_data.items())
    print(f"  Original sequences: {total_original}")
    print(f"  De-duplicated sequences: {len(mapping_data)}")
    print(f"  Reduction: {(1 - len(mapping_data)/total_original)*100:.1f}%")

def main():
    args = parse_args()
    
    # Read genome labels if provided
    genome_labels = {}
    if args.genome_labels:
        genome_labels = read_genome_labels(args.genome_labels)
    
    # Perform genome-aware deduplication
    dedup_fasta_genome_aware(
        args.input, 
        args.output, 
        args.mapping,
        genome_labels, 
        args.separator
    )

if __name__ == '__main__':
    main()
