#!/usr/local/python

import gzip
import json
import subprocess

output_fp = "/usr/map_viruses/tests/example.fastq.json.gz"

subprocess.call([
    "map_viruses.py",
    "--input", "/usr/map_viruses/tests/example.fastq",
    "--metadata", "/usr/map_viruses/tests/example.tsv",
    "--ref-db", "/usr/map_viruses/tests/example.dmnd",
    "--output-path", "/usr/map_viruses/tests/example.fastq.json.gz",
])

output = json.load(gzip.open(output_fp))

protein_abund = output["results"]["proteins"]
genome_dat = output["results"]["genomes"]

assert len(protein_abund) == 11

assert len(genome_dat) == 1

assert genome_dat[0]["genome"] == "NC_001422.1"
assert genome_dat[0]["nreads"] == 20833
assert genome_dat[0]["total_proteins"] == 11
assert genome_dat[0]["detected_proteins"] == 11
assert genome_dat[0]["total_length"] == 2327

print("Success")
