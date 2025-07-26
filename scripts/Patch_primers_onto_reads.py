#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Mapping for IUPAC ambiguity codes to lexically smallest base
iupac_lexical_min = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'A', 'Y': 'C', 'S': 'C', 'W': 'A',
    'K': 'G', 'M': 'A', 'B': 'C', 'D': 'A',
    'H': 'A', 'V': 'A', 'N': 'A'
}

def get_lexical_seq(seq):
    return ''.join(iupac_lexical_min.get(base.upper(), 'A') for base in seq)

def main():
    parser = argparse.ArgumentParser(description="Patch synthetic primers onto reads")
    parser.add_argument('--reads', required=True, help='Input FASTA file of reads')
    parser.add_argument('--primers', required=True, help='2-line text file: forward, reverse primer')
    parser.add_argument('--output', required=True, help='Output patched FASTA file')
    args = parser.parse_args()

    # Read primer sequences
    with open(args.primers) as pf:
        lines = [line.strip() for line in pf if line.strip()]
        if len(lines) != 2:
            raise ValueError("Primers file must contain exactly 2 non-empty lines")
        fwd_primer = get_lexical_seq(lines[0])
        rev_primer = get_lexical_seq(lines[1])
        rev_primer_rc = str(Seq(rev_primer).reverse_complement())

    # Patch each sequence
    with open(args.output, 'w') as out_handle:
        for record in SeqIO.parse(args.reads, "fasta"):
            patched_seq = fwd_primer + str(record.seq) + rev_primer_rc
            record.seq = Seq(patched_seq)
            SeqIO.write(record, out_handle, "fasta")

if __name__ == '__main__':
    main()

