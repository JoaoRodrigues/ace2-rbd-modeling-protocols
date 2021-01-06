#!/usr/bin/env python

import re

from Bio import SeqIO


species_re = re.compile("\[([\w\s]+)\]")

# Read species list
with open("species.list") as handle:
    species_set = {l.strip() for l in handle if l[0] != "#"}
print(f"Read {len(species_set)} species from selection file.")

# Build seq database
seqdict = {}  # sp name -> seq

with open("ace2-orthologs-all.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.description in seqdict:
            print(f"Duplicated sequence: {record.id}")
            continue

        record_species = species_re.findall(record.description)
        if record_species and record_species[0] in species_set:
            species_name = "_".join(record_species[0].strip().split())
            seqdict[species_name] = str(record.seq)

with open("ace2-orthologs-extra.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.description in seqdict:
            print(f"Duplicated sequence: {record.id}")
            continue

        record_species = species_re.findall(record.description)
        if record_species and record_species[0] in species_set:
            species_name = "_".join(record_species[0].strip().split())
            seqdict[species_name] = str(record.seq)

print(f"Selected {len(seqdict)} sequences")

# Filter dataset
with open("ace2-orthologs-selected.fasta", "w") as handle:
    for name, seq in seqdict.items():
        print(f">{name}\n{seq}", file=handle)
