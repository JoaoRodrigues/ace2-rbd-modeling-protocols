#!/usr/bin/env python

import argparse
from Bio import AlignIO


ap = argparse.ArgumentParser()
ap.add_argument("input")
ap.add_argument("output")
args = ap.parse_args()

# Read
with open(args.input) as handle:

    aln = AlignIO.read(handle, "fasta")

    # Find indices in gapped reference
    refseq = aln[0]
    assert refseq.description == "6M0J", "Wrong reference."

    boundaries = []
    for idx, aa in enumerate(refseq):
        if aa != '-':
            boundaries.append(idx)

    # Now trim alignment
    sta, end = boundaries[0], boundaries[-1] + 1
    seqdict = {}
    for record in aln:
        seq = str(record.seq)[sta:end]
        name = record.description
        seqdict[name] = seq


# Write
with open(args.output, "w") as handle:
    for name, seq in seqdict.items():
        print(f">{name}\n{seq}", file=handle)
