#!/usr/bin/env python
"""Pre-processes the initial sequence file.

* Shortens the description for each sequence.
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
GISAID_LABELS = [
    "gene_name", "isolate_name", "date", "isolate_id", "passage_details",
    "type_location", "host", "origin_lab", "submitter_lab", "submitter",
    "country"
]
# Extra labels not in the GISAID file.
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
        help="Input file in FASTA format."
    )
    ap.add_argument(
        "output",
        type=Path,
        help="Path to output file."
    )
    ap.add_argument(
        "--min-length",
        type=int,
        default=100,
        help="Minimum sequence length to consider."
    )
    ap.add_argument(
        "--include-ambig",
        action="store_true",
        help="Include ambiguous characters in length calculation."
    )

    return ap.parse_args()


def read_seqfile(handle, min_length=100, include_ambig=False):
    """Read sequences in FASTA format into two dictionaries.
    
    One contains the sequences, another metadata.
    """

    # Define length function
    # Returns True if the sequence *passes* the length check.
    if include_ambig:
        check_length = lambda s: len(s) > min_length
    else:
        check_length = lambda s: (len(s) - s.count("X")) > min_length

    # Store metadata in a list of dicts
    metadata = []

    seqdict = {}
    parsedset = set()  # to catch duplicated isolate ids

    seqiter = SeqIO.parse(handle, "fasta")

    n_short = 0
    for idx, record in enumerate(seqiter, start=1):
        seq = str(record.seq)

        tokens = record.description.split("|")
        assert len(tokens) == len(GISAID_LABELS), \
            "Entry on line {idx} does have the expected description fields."

        datadict = dict(zip(GISAID_LABELS, tokens))

        name = datadict["isolate_id"] + "|" + datadict["host"]

        if name in parsedset:
            logging.warning(
                f"Duplicated entry '{name}': skipping."
            )
            continue
        
        parsedset.add(name)

        # Add basic sequence statistics
        datadict["seq_length"] = len(seq)
        datadict["num_ambig"] = seq.count("X")

        is_short = not check_length(seq)

        # Catch redundant sequences
        in_seqdict = seqdict.get(seq)
        if in_seqdict is None:
            datadict["representative"] = name  # itself
            if not is_short:
                datadict["is_filtered"] = False
                seqdict[seq] = name
            else:
                datadict["is_filtered"] = True
                n_short += 1
        else:
            datadict["representative"] = in_seqdict
            datadict["is_filtered"] = False
        
        metadata.append(datadict)

    logging.info(
        f"Read {len(seqdict)} unique sequences out of {len(metadata)} total."
    )
    logging.info(
        f"Discarded {n_short} short sequences."
    )

    return (seqdict, metadata)


if __name__ == "__main__":

    args = read_args()

    # Read sequence data
    with open(args.aln) as handle:
        seqdict, metadata = read_seqfile(
            handle,
            min_length=args.min_length,
            include_ambig=args.include_ambig
        )

    # Write non-redundant sequences as FASTA
    logging.info(f"Writing non-redundant sequences to {args.output}")
    with open(args.output, "w") as handle:
        for seq, name in seqdict.items():
            print(f">{name}\n{seq}", file=handle)

    # Write metadata as csv file.
    logging.info(f"Writing sequence metadata to metadata.csv")
    with open("metadata.csv", "w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=GISAID_LABELS+EXTRA_LABELS,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(metadata)

    logging.info("Done.")