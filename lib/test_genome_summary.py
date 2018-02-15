#!/usr/bin/python

import pandas as pd
from aln_helpers import parse_alignment
from aln_helpers import summarize_genomes

fp = "/usr/map_viruses/tests/example.aln"
metadata_fp = "/usr/map_viruses/tests/example.tsv"

protein_abund = parse_alignment(fp)

metadata = pd.read_table(metadata_fp, sep='\t')

protein_abund, genome_dat = summarize_genomes(protein_abund, metadata)

assert len(protein_abund) == 11

assert len(genome_dat) == 1

assert genome_dat[0]["genome"] == "NC_001422.1"
assert genome_dat[0]["nreads"] == 20833
assert genome_dat[0]["total_proteins"] == 11
assert genome_dat[0]["detected_proteins"] == 11
assert genome_dat[0]["total_length"] == 2327

print("Success")
