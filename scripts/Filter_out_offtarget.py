#!/usr/bin/env python3

import sys

def filter_fasta(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        write_seq = False
        header = ''
        seq_lines = []
        
        for line in infile:
            line = line.rstrip()
            if line.startswith('>'):
                # Write previous sequence if flagged
                if write_seq and header:
                    outfile.write(f"{header}\n{''.join(seq_lines)}\n")
                # Check orientation
                if '.orientation="FF"' in line or '.orientation="RR"' in line:
                    write_seq = False
                else:
                    write_seq = True
                    header = line
                    seq_lines = []
            else:
                if write_seq:
                    seq_lines.append(line)

        # Write last sequence if needed
        if write_seq and header:
            outfile.write(f"{header}\n{''.join(seq_lines)}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python filter_fasta.py input.fasta output.fasta")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    filter_fasta(input_fasta, output_fasta)

