# Processing raw GISAID Spike protein dataset.
This folder contains scripts to align and process raw sequence data from GISAID.

## How to run
Run `make_aln.sh <seqfile>` where seqfile is the spike protein data downloaded
from GISAID.

## Rationale
### Alignment
The script will remove short sequences and create a non-redundant (100% seq. id)
set of sequences that will be aligned by MAFFT. The remaining non-redundant sequences
can be added back to the alignment at the end. Finally, the alignment is "cleaned",
by removing sequences with low coverage (gaps) and trimming to the region defined in
the crystal structure (6M0J).

### Post-Processing
The non-redundant alignment will be filtered for coverage and also processed for
ambiguous characters. Since there is no good way of handling Xs in these sequences we
choose to revert them to the amino acids in the reference sequence. This yields ~100
more variants than simply discarding these sequences and allows us to model them.

### Re-alignment
The set of non-redundant RBD-only sequences are then re-aligned using MAFFT and a slow
iterative refinement alignment method (G-INS-i) to fix any issues of with the previous
alignment due to ambig characters.

