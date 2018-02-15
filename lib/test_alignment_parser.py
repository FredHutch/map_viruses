#!/usr/bin/python

from aln_helpers import parse_alignment

fp = "/usr/map_viruses/tests/example.aln"

protein_abund = parse_alignment(fp)

# There are 11 proteins present
assert len(protein_abund) == 11

for prot in protein_abund:
    for k in ["coverage", "bitscore", "depth", "pctid", "alen"]:
        assert prot[k] > 0
        assert isinstance(prot[k], float)
    for k in ["nreads", "length"]:
        assert prot[k] > 0
        assert isinstance(prot[k], int)
    assert isinstance(prot["protein"], str)

print("Success")
