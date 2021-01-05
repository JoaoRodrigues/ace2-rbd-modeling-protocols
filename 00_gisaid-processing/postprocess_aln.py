#!/usr/bin/env python
"""Processes aligned sequences.

* Trims alignment to a specific region (reference sequence)
* Removes short sequences.
* Removes redundancy.
* Backmaps ambiguous characters to match reference sequence.
* Produces a CSV file with metadata for each entry.
"""

import argparse
import csv
import logging
from pathlib import Path
import sys

from Bio import SeqIO


# labels from GISAID readme (as of Dec 2020):
# >Gene name|Isolate name|YYYY-MM-DD|Isolate ID|Passage details/history|
# Type^^location/state|Host|Originating lab|Submitting lab|Submitter|
# Location(country)
LABELS = ["isolate_id", "host"]
# Extra labels added during processing
EXTRA_LABELS = ["seq_length", "num_ambig", "representative", "is_filtered"]


logging.basicConfig(
    format="[%(asctime)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)


def read_args():
    """Read command-line arguments."""

    logging.info(f"cmd: {' '.join(sys.argv)}")

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "aln",
        type=Path,
        help="Input alignment in FASTA format."
    )
    ap.add_argument(
        "output",
        type=Path,
        help="Path to output file."
    )
    # Filters
    ap.add_argument(
        "--max-ambig-fraction",
        type=float,
        default=0.5,
        help="Maximum ambiguous characters as a fraction of sequence length."
    )
    ap.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        help="Minimum coverage as a fraction of reference sequence length."
    )

    return ap.parse_args()


def read_aln(handle):
    """Read and trim aligned sequences."""

    # Read alignment into dictionary
    with open(args.aln) as handle:
        seqdict = {
            r.description: str(r.seq)
            for r in SeqIO.parse(handle, "fasta")
        }
        logging.info(f"Read {len(seqdict)} aligned sequences.")

    refname = next(iter(seqdict))
    logging.info(f"Using as reference sequence: {refname}")
    refseq = seqdict[refname]

    # Boundaries to trim alignment
    idxsort = sorted(
        idx for idx, aa in enumerate(refseq) if aa != "-"
    )
    beg, end = idxsort[0], idxsort[-1] + 1
    logging.info(f"Trimming sequences to region {beg}-{end - 1}")

    # Trim alignment
    for name in seqdict:
        seqdict[name] = seqdict[name][beg:end]

    return seqdict


def process_aligned_sequences(aln, frac_max_ambig=0.5, min_coverage=0.8):
    """Process the alignment to remove redundancy and low-quality entries."""

    alndict = {}  # store postprocessed data.
    metadata = [] # store metadata in a list of dicts

    refname = next(iter(aln))
    refseq = aln[refname]

    ref_nongaps = {idx for idx, aa in enumerate(refseq) if aa != "-"}

    n_skip = 0
    for name, seq in aln.items():
        skip = False

        tokens = name.split("|")
        assert len(tokens) == len(LABELS), \
            "Entry {name} does have the expected description fields."

        datadict = dict(zip(LABELS, tokens))
        datadict["is_filtered"] = False

        # Add basic sequence statistics
        datadict["seq_length"] = len(seq)
        datadict["num_ambig"] = seq.count("X")
        datadict["num_gaps"] = seq.count("-")

        # Discard sequences with more than X ambiguous characters
        if seq.count("X") / len(seq) > frac_max_ambig:
            datadict["is_filtered"] = True
            skip = True

        # Discard sequences with low coverage of the reference
        # Coverage is calculated as the fraction of ungapped positions
        # matching ungapped reference sites.
        n_refmatched_gaps = sum(
            1 for idx, aa in enumerate(seq)
            if idx in ref_nongaps and aa == "-"
        )
        if 1 - n_refmatched_gaps / len(ref_nongaps) < 0.8:
            datadict["is_filtered"] = True
            skip = True

        # Backmap to reference
        # This has the side effect of backmapping Xs to gaps as well,
        # which might mess up the alignment if there are insertions.
        # Thus, we re-align more carefully in a later step.
        if not skip:
            seq = "".join(
                aa if aa != "X" else refseq[idx]
                for idx, aa in enumerate(seq)
            )

        # Catch redundant sequences
        # We might have gaps at termini, which can be a sequencing
        # artifact. Let's trim those before checking for redundancy.
        in_dict = alndict.get(seq.strip("-"))

        if in_dict is None:
            datadict["representative"] = name  # itself
            if not skip:
                alndict[seq] = name
            else:
                n_skip += 1
        else:
            datadict["representative"] = in_dict
        
        metadata.append(datadict)

    logging.info(
        f"Read {len(alndict)} unique sequences out of {len(metadata)} total."
    )
    logging.info(
        f"Discarded {n_skip} low-quality sequences."
    )

    # Invert dict - safe bc we have unique isolate_ids.
    alndict = {val: key for key, val in alndict.items()}

    # Finally, remove fully gapped columns
    trseq = zip(*alndict.values())
    nseqs = len(alndict)
    rmcol = []
    for idx, col in enumerate(trseq):
        if col.count("-") == nseqs:
            rmcol.append(idx)

    if rmcol:
        logging.info(f"Will remove {len(rmcol)} 100%-gapped columns.")

        rmcol = set(rmcol)
        for name in list(alndict):
            alndict[name] = "".join(
                aa for aaidx, aa in enumerate(alndict[name])
                if aaidx not in rmcol
            )

    return (alndict, metadata)


if __name__ == "__main__":

    args = read_args()

    # Read sequence data
    with open(args.aln) as handle:
        aln = read_aln(handle)

    alndict, metadata = process_aligned_sequences(
        aln,
        args.max_ambig_fraction,
        args.min_coverage
    )

    # Write non-redundant sequences as FASTA
    logging.info(f"Writing non-redundant sequences to {args.output}")
    with open(args.output, "w") as handle:
        for name, seq in alndict.items():
            print(f">{name}\n{seq}", file=handle)

    # Write metadata as csv file.
    logging.info(f"Writing alignment metadata to metadata_aln.csv")
    with open("metadata_aln.csv", "w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=LABELS+EXTRA_LABELS,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(metadata)

    logging.info("Done.")