#!/usr/bin/python
"""Make a database of viral genomes."""

import gzip
import logging
import argparse
import pandas as pd
from Bio import GenBank
from lib.exec_helpers import run_cmds


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Download and format a reference database of viral genomes.
    """)

    parser.add_argument("--prefix",
                        type=str,
                        help="""Prefix for output files.""")

    args = parser.parse_args()

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [map_viruses.py] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write logs to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    dat = []

    # Download viral genomes from NCBI
    for url in [
        "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.gpff.gz",
        "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.gpff.gz",
    ]:
        run_cmds([
            "wget", url
        ])
        fp = url.split("/")[-1]
        with gzip.open(fp, "rt") as handle:
            for ix, record in enumerate(GenBank.parse(handle)):
                taxid = None
                product = None
                locus_tag = None
                coded_by = None
                genome = None
                genome_range = None

                for feature in record.features:
                    for qualifier in feature.qualifiers:
                        if feature.key == "source":
                            if qualifier.key == "/db_xref=":
                                taxid = qualifier.value.strip('"')
                                taxid = taxid.replace("taxon:", "")
                        elif feature.key == "Protein":
                            if qualifier.key == "/product=":
                                product = qualifier.value.strip('"')
                        elif feature.key == "CDS":
                            if qualifier.key == "/locus_tag=":
                                locus_tag = qualifier.value.strip('"')
                            elif qualifier.key == "/coded_by=":
                                coded_by = qualifier.value.strip('"')
                                genome, genome_range = coded_by.split(":", 1)
                                genome = genome.replace("complement(", "")
                                genome = genome.replace("join(", "")
                dat.append({
                    "protein": record.locus,
                    "sequence": record.sequence,
                    "organism": record.organism,
                    "taxonomy": '; '.join(record.taxonomy),
                    "definition": record.definition,
                    "taxid": taxid,
                    "product": product,
                    "locus_tag": locus_tag,
                    "genome": genome,
                    "region": genome_range,
                    "length": len(record.sequence),
                })
                if ix > 0 and ix % 1e3 == 0:
                    logging.info("Parsed {:,} sequence records".format(ix))

    df = pd.DataFrame(dat)

    # Make sure all of the protein names are unique
    assert df.shape[0] == len(df["protein"].unique())

    logging.info("Writing protein sequences to {}.fastp".format(args.prefix))
    with open("{}.fastp".format(args.prefix), "wt") as handle:
        for ix, r in df.iterrows():
            handle.write(">{}\n{}\n".format(r["protein"], r["sequence"]))

    del df["sequence"]

    logging.info("Writing mappings to {}.tsv".format(args.prefix))
    df.to_csv("{}.tsv".format(args.prefix), sep="\t", index=None)

    logging.info("Formatting the DIAMOND database")
    run_cmds([
        "diamond", "makedb",
        "--in", "{}.fastp".format(args.prefix),
        "--db", args.prefix
    ])
