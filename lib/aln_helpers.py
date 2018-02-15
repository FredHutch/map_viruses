#!/usr/bin/python

import logging
import numpy as np
from collections import defaultdict


def read_mapping_tsv(mapping_fp):
    """Read in the mapping file."""
    logging.info("Reading in the mapping file ({})".format(mapping_fp))
    with open(mapping_fp, "rt") as f:
        header = f.readline()
    header = header.rstrip("\n").split("\t")
    assert 'protein' in header
    assert 'genome' in header
    assert 'length' in header

    # Keep track of all of the information for a given genome
    genome_dat = defaultdict(lambda: defaultdict(set))

    # Keep track of which protein is a part of which genome
    protein_genome = {}

    # Keep track of the length of each protein
    protein_length = {}

    with open(mapping_fp, "rt") as f:
        for ix, line in enumerate(f):
            if len(line) == 1:
                continue
            elif ix == 0:
                continue
            else:
                line = line.rstrip("\n").split("\t")
                assert len(line) == len(header)

                entry = dict(zip(header, line))

                protein_genome[entry["protein"]] = entry["genome"]
                protein_length[entry["protein"]] = entry["length"]

                for k, v in entry.items():
                    if k not in ["protein", "genome", "length"]:
                        genome_dat[entry["genome"]][k].add(v)

    # Reformat the genome dat as a dict of dicts, with values as strings
    genome_dat = {
        genome: {
            k: "; ".join(list(v))
            for k, v in values.items()
        }
        for genome, values in genome_dat.items()
    }

    return protein_genome, protein_length, genome_dat


def summarize_genomes(protein_abund, mapping):
    """From a set of protein abundances, summarize the genomes."""
    protein_genome, protein_length, genome_dat = mapping

    return


def parse_alignment(align_fp,
                    subject_len,
                    subject_ix=1,
                    pctid_ix=2,
                    alen_ix=3,
                    sstart_ix=6,
                    send_ix=7,
                    bitscore_ix=9):
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
                assert s in subject_len
                coverage[s] = np.zeros(subject_len[s], dtype=int)

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
