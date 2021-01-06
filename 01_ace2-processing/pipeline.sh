#!/usr/bin/env bash

set -e

rm -f ace2-orthologs-selected.*
rm -f ace2.final.fasta

python make_selection.py
./make_aln.sh ace2-orthologs-selected.fasta

# Trim to template
python trim_aln.py ace2-orthologs-selected.aln.fasta ace2-orthologs-selected.aln.final.fasta
