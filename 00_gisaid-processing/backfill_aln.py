#!/usr/bin/env python
"""Adds unaligned redundant sequences to a non-redundant alignment."""

import argparse
import collections
from pathlib import Path
import sys

from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument(
    "alnfile",
    type=Path,
    help="Input alignment file in FASTA format."
)
ap.add_argument(
    "seqfile",
    type=Path,
    help="Input sequence file in FASTA format."
)
ap.add_argument(
    "outfile",
    type=Path,
    help="Output sequence file name."
)

args = ap.parse_args()

print(f"cmd: {' '.join(sys.argv)}")

# Read and map sequence to ids
with open(args.seqfile) as handle:

    seqdict = collections.defaultdict(list)
    
    n_read = 0
    for record in SeqIO.parse(handle, "fasta"):
        n_read += 1

        seq = str(record.seq)
        name = record.description
        seqdict[seq].append(name)

print(f"Read {n_read} sequences ({len(seqdict)} unique)")

# Read in alignment
with open(args.alnfile) as handle:
    alndict = {}

    for record in SeqIO.parse(handle, "fasta"):
        seq = str(record.seq)
        name = record.description
        alndict[name] = seq

print(f"Read {len(alndict)} aligned sequences.")

# Now match names in seqdict to names in alndict
seqgroups = list(seqdict.values())
matchdict = {}
for n in alndict.keys():
    for g in seqgroups:
        if n in set(g):
            break
    else:
        print(f"Could not match {n} to original seqfile.")
        print("This is OK if this is the reference.")
        g = [n]

    matchdict[n] = g
    try:
        seqgroups.remove(g)
    except:
        pass  # reference is not in seqgroups

# Write aligned sequences in order of seqfile.
with open(args.outfile, "w") as handle:
    for alnkey in matchdict:
        namelist = matchdict[alnkey]
        seq = alndict[alnkey]
        
        for name in namelist:
            print(f">{name}\n{seq}", file=handle)
print("Done")