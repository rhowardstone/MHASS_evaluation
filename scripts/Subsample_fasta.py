#!/usr/bin/env python3
"""
Subsample_fasta.py
Randomly subsample up to N sequences from a (possibly huge) FASTA/FASTA.GZ file,
uniformly without replacement, using reservoir sampling.

CLI (compatible with existing pipeline call):
  -i / --input   : input FASTA(.gz)
  -o / --output  : output FASTA
  -n / --num_seqs: max sequences to keep (default 100)
  -r / --random  : kept for backward compatibility; ignored (sampling is always random)
  -s / --seed    : PRNG seed (default 42)
"""

import argparse
import gzip
import random
import sys
from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Uniformly subsample sequences from a FASTA file "
                    "without loading the whole file into memory."
    )
    p.add_argument("-i", "--input", required=True, help="Input FASTA or FASTA.GZ")
    p.add_argument("-o", "--output", required=True, help="Output FASTA")
    p.add_argument("-n", "--num_seqs", type=int, default=100,
                   help="Maximum number of sequences to sample (default: 100)")
    # kept for pipeline compatibility
    p.add_argument("-r", "--random", action="store_true",
                   help="(Ignored — sampling is always random.)")
    p.add_argument("-s", "--seed", type=int, default=42,
                   help="Random seed (default: 42)")
    return p.parse_args()


def open_maybe_gzip(path: Path, mode: str = "rt"):
    """Transparent open for plain or gzipped files."""
    return gzip.open(path, mode) if str(path).endswith((".gz", ".gzip")) else open(path, mode)


def reservoir_sample(seq_iter, k: int) -> List[Tuple[int, SeqRecord]]:
    """
    Reservoir-sample k items from seq_iter.
    Returns list of (original_index, SeqRecord) so we can restore input order.
    """
    reservoir = []  # type: List[Tuple[int, SeqRecord]]
    for idx, record in enumerate(seq_iter):
        if idx < k:
            reservoir.append((idx, record))
            continue
        j = random.randint(0, idx)
        if j < k:
            reservoir[j] = (idx, record)
    return reservoir


def main() -> None:
    args = parse_args()
    random.seed(args.seed)

    in_path = Path(args.input)
    out_path = Path(args.output)

    try:
        with open_maybe_gzip(in_path, "rt") as handle:
            # stream through the file once
            sampled_pairs = reservoir_sample(SeqIO.parse(handle, "fasta"), args.num_seqs)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to read '{in_path}': {e}\n")
        sys.exit(1)

    if not sampled_pairs:  # empty file
        sys.stderr.write(f"[WARNING] No sequences found in '{in_path}'. Output left empty.\n")
        open(out_path, "w").close()
        return

    # Restore original order before writing
    sampled_pairs.sort(key=lambda p: p[0])
    sampled_records = [rec for _, rec in sampled_pairs]

    try:
        with open(out_path, "w") as out_handle:
            SeqIO.write(sampled_records, out_handle, "fasta")
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to write '{out_path}': {e}\n")
        sys.exit(1)

    sys.stderr.write(
        f"Subsampled {len(sampled_records)} / {max(sampled_pairs)[0] + 1} "
        f"sequences from '{in_path.name}' → '{out_path.name}'\n"
    )


if __name__ == "__main__":
    main()

