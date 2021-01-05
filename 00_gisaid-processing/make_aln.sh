#!/usr/bin/env bash

# Align Spike protein sequences from GISAID sequences.
# Alignment in two steps:
#   First, fast align all full-length non-redundant sequences
#   Then, trim to region of template/reference, backmap, and slow(er) re-align.
#

MAFFT="/home/joaor/software/mafft_nov2020/bin/mafft"

seqfile=$1
if [ ! -f "$seqfile" ]
then
  echo "usage: $0 <sequence file>"
  echo "example: $0 spikeprot0101.fasta"
  exit 1
else
  rootname=${seqfile%%.fasta}
fi

rm -f ${rootname}.aln.*
rm -f ${rootname}.prealn.*

# Pre-process sequences before alignment. Remove redundancy to
# allow faster alignment. We can backfill later if necessary.
# As of Dec 2020, reduces alignment from 250k to ~40k sequences.
python preprocess_sequences.py $seqfile ${rootname}.preproc.fasta

# Add reference RBD-only sequence on top for mapping later.
cat SPIKE_6M0J.fasta ${rootname}.preproc.fasta > ${rootname}.prealn.fasta

# Align non-redundant sequences using automatic MAFFT (FFT-NS-2 most likely)
# Use BLOSUM80 (closely related viral sequences)
mafft_cmd="$MAFFT --auto --thread 24 --bl 80"
echo "cmd: $mafft_cmd"
$mafft_cmd ${rootname}.prealn.fasta > ${rootname}.aln.fasta

# Post-process the alignment
#   1. Trim to reference region
#   2. Remove short sequences (coverage/ambig).
#   2. Remove redundant sequences
#   3. Backmap ambiguous characters.
python postprocess_aln.py ${rootname}.aln.fasta ${rootname}.aln.postproc.fasta

# Finally, re-align carefully with an accurate method.
# This avoids issues with prev ambig characters.
mafft_cmd="$MAFFT --thread 24 --bl 80 --globalpair --maxiterate 16 --inputorder"
echo "cmd: $mafft_cmd"
$mafft_cmd ${rootname}.aln.postproc.fasta > ${rootname}.aln.final.fasta

# Uncomment to remove temporary files.
#rm -f ${rootname}.aln.preproc.fasta
#rm -f ${rootname}.prealn.fasta
#rm -f ${rootname}.aln.postproc.fasta