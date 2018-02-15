#!/usr/bin/python

import logging
import numpy as np
from collections import defaultdict


def summarize_genomes(protein_abund, mapping):
    """From a set of protein abundances, summarize the genomes."""

    # Cache the stats for the proteins
    prot_stats = {
        k: {
            p["protein"]: p[k]
            for p in protein_abund
        }
        for k in ["coverage", "depth", "pctid", "bitscore", "alen", "nreads"]
    }

    genome_abund = []
    # Iterate over all of the genomes
    for genome, proteins in mapping.groupby("genome"):

        # Calculate the aggregate coverage, depth, number of proteins, etc.

        # Number of proteins with any reads detected
        cov_prots = sum(proteins["protein"].apply(lambda p: p in prot_stats))
        if cov_prots == 0:
            continue
        logging.info("Collecting stats for {}".format(genome))
        # Total number of proteins in this genome
        tot_prots = proteins.shape[0]

        dat = {
            "genome": genome,
            "total_proteins": tot_prots,
            "detected_proteins": cov_prots,
        }

        # Add stats to the table
        for stat in [
            "coverage", "depth", "pctid", "bitscore", "alen", "nreads"
        ]:
            proteins[stat] = proteins["protein"].apply(
                lambda p: prot_stats[stat].get(p, 0))
            dat[stat] = (proteins[stat] * proteins["length"]).sum() / proteins["length"].sum()

        genome_abund.append(dat)

    return genome_abund


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
    pctid = defaultdict(list)
    alen = defaultdict(list)
    bitscore = defaultdict(list)

    logging.info("Reading from {}".format(align_fp))
    with open(align_fp, "rt") as f:
        for ix, line in enumerate(f):
            line = line.rstrip("\n").split("\t")
            s = line[subject_ix]
            if s not in coverage:
                slen = int(line[slen_ix])
                coverage[s] = np.zeros(slen, dtype=int)

            pctid[s].append(float(line[pctid_ix]))
            alen[s].append(int(line[alen_ix]))
            bitscore[s].append(float(line[bitscore_ix]))

            coverage[s][
                int(line[sstart_ix]): int(line[send_ix])
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
            "nreads": len(pctid[s])
        })
        if ix > 0 and ix % 1e3 == 0:
            logging.info("Summarized coverage for {:,} subjects".format(ix))

    logging.info("Summarized coverage for {:,} subjects".format(ix))

    return output
