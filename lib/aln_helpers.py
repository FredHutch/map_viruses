#!/usr/bin/python

import logging
import numpy as np
import pandas as pd
from collections import defaultdict


def summarize_genomes(protein_abund, metadata):
    """From a set of protein abundances, summarize the genomes."""

    # Format the protein data as a DataFrame
    protein_abund = pd.DataFrame(protein_abund).set_index("protein")

    # Add the detected protein information to the metadata table
    for k in [
        "coverage", "depth", "pctid", "bitscore", "alen", "nreads"
    ]:
        metadata[k] = metadata["protein"].apply(
            protein_abund[k].to_dict().get
        ).fillna(0)

    assert (metadata["length"] > 0).all()

    # Subset to the GENOMES that have _any_ proteins detected
    metadata = pd.concat([
        genome_dat
        for genome, genome_dat in metadata.groupby("genome")
        if (genome_dat["coverage"] > 0).any()
    ])

    # Save the protein summary for these genomes
    protein_abund = metadata.to_dict(orient="records")

    # Now make a summary on a per-genome basis
    genome_abund = []
    # Iterate over all of the genomes
    for genome, proteins in metadata.groupby("genome"):

        # Calculate the aggregate genome length
        agg_len = proteins["length"].sum()

        # Calculate the aggregate coverage, depth, number of proteins, etc.
        dat = {
            "total_length": int(agg_len),
            "total_proteins": proteins.shape[0],
            "detected_proteins": (proteins["coverage"] > 0).sum(),
            "genome": genome,
            "nreads": int(proteins["nreads"].sum()),
        }

        for k in ["coverage", "depth", "pctid", "bitscore", "alen"]:
            # Make a length-adjusted average
            dat[k] = (proteins[k] * proteins["length"]).sum() / agg_len

        # For all of the other columns, add them if they are unique
        for k in proteins.columns:
            if k not in dat:
                if len(proteins[k].unique()) == 0:
                    dat[k] = proteins[k].values[0]

        genome_abund.append(dat)

    return protein_abund, genome_abund


def parse_alignment(align_fp,
                    subject_ix=1,
                    pctid_ix=2,
                    alen_ix=3,
                    sstart_ix=6,
                    send_ix=7,
                    bitscore_ix=9,
                    slen_ix=11):
    """
    Parse an alignment in BLAST6 format and calculate coverage per subject.
    """

    # Keep track of a number of different metrics for each subject
    coverage = {}
    subject_len = {}
    pctid = defaultdict(list)
    alen = defaultdict(list)
    bitscore = defaultdict(list)

    logging.info("Reading from {}".format(align_fp))
    with open(align_fp, "rt") as f:
        for ix, line in enumerate(f):
            if len(line) == 1:
                continue
            line = line.rstrip("\n").split("\t")
            s = line[subject_ix]
            if s not in coverage:
                slen = int(line[slen_ix])
                subject_len[s] = slen
                coverage[s] = np.zeros(slen, dtype=int)

            pctid[s].append(float(line[pctid_ix]))
            alen[s].append(int(line[alen_ix]))
            bitscore[s].append(float(line[bitscore_ix]))

            coverage[s][
                (int(line[sstart_ix]) - 1): int(line[send_ix])
            ] += 1

            if ix > 0 and ix % 1e6 == 0:
                logging.info("Parsed {:,} alignments".format(ix))

    logging.info("Parsed {:,} alignments".format(ix))

    # Calculate the per-subject stats
    output = []
    for ix, s in enumerate(coverage):
        output.append({
            "protein": s,
            "coverage": (coverage[s] > 0).mean(),
            "depth": coverage[s].mean(),
            "pctid": np.mean(pctid[s]),
            "alen": np.mean(alen[s]),
            "bitscore": np.mean(bitscore[s]),
            "nreads": len(pctid[s]),
            "length": subject_len[s],
        })
        if ix > 0 and ix % 1e3 == 0:
            logging.info("Summarized coverage for {:,} subjects".format(ix))

    logging.info("Summarized coverage for {:,} subjects".format(ix))

    return output
