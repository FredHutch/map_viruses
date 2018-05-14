#!/usr/bin/env python

import os
import gzip
import argparse
import pandas as pd
from subprocess import call
from Bio.SeqIO.FastaIO import SimpleFastaParser

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Reformat the IMG/VR reference data for use with the `map_viruses` tool.
    """)

    parser.add_argument("--img-vr-metadata",
                        type=str,
                        required=True,
                        help="""Metadata file provided by IMG/VR.""")
    parser.add_argument("--img-vr-proteins",
                        type=str,
                        required=True,
                        help="""Protein sequence (faa) file provided by IMG/VR.""")
    parser.add_argument("--output-prefix",
                        type=str,
                        required=True,
                        help="Prefix for output files.")

    args = parser.parse_args()

    assert os.path.exists(args.img_vr_metadata)
    assert os.path.exists(args.img_vr_proteins)

    # Read in the metadata
    print("Reading in " + args.img_vr_metadata)
    metadata = pd.read_table(args.img_vr_metadata)
    assert metadata.shape[0] == metadata["mVCs"].unique().shape[0]

    # Keep track of all of the genome IDs with metadata
    all_genome_ids = set(metadata["mVCs"].tolist())

    # Index on the genome ID
    metadata.set_index("mVCs", inplace=True)

    # Read in each of the protein names and genome names from the FAA input file
    proteins = []
    print("Reading in " + args.img_vr_proteins)
    with gzip.open(args.img_vr_proteins, "rt") as f:
        ix = 0
        for header, seq in SimpleFastaParser(f):
            # Split the header by "_____"
            genome_id = "_____".join(header.split("_____")[:2])
            # Make sure we have metadata for this genome
            assert genome_id in all_genome_ids, "No metadata for " + genome_id

            # Add to the list
            proteins.append({
                "protein": header,
                "length": len(seq),
                "genome": genome_id
            })
            ix += 1
            if ix % 100000 == 0:
                print("Read in {:,} proteins".format(ix))

    # Format as a DataFrame
    proteins = pd.DataFrame(proteins)

    # Now add the metadata from IMG/VR
    for k in metadata.columns:
        print("Adding metadata for " + k)
        proteins[k] = proteins["genome"].apply(metadata[k].to_dict().get)

    proteins.to_csv(args.output_prefix + ".tsv", sep="\t", index=False)

    # Now index the protein sequences with DIAMOND
    db_fp = args.output_prefix + ".dmnd"
    print("Making DIAMOND database, writing to " + db_fp)
    call(["diamond", "makedb", "--in", args.img_vr_proteins, "--db", db_fp])
    assert os.path.exists(db_fp)
