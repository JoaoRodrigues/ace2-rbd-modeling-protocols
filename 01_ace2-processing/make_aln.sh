#!/usr/bin/env bash

MAFFT="/home/joaor/software/mafft_nov2020/bin/mafft"


seqfile=$1
if [ ! -f "$seqfile" ]
then
  echo "usage: $0 <sequence file>"
  exit 1
else
  rootname=${seqfile%%.fasta}
fi

rm -f ${rootname}.aln.*
rm -f ${rootname}.prealn.*

cat ACE2_6M0J.fasta $1 > ${rootname}.prealn.fasta

# Accurate alignment
$MAFFT --thread 24 --globalpair --maxiterate 16 --inputorder ${rootname}.prealn.fasta > ${rootname}.aln.fasta
