# Map Viruses
### Align short reads against viral genomes in protein space

[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/map_viruses/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/map_viruses)


Researchers frequently want to detect viral genomes from short read metagenomic datasets. 
My preferred way of doing this is to align reads against those viral genomes in amino acid
space (akin to BLASTx) using the DIAMOND aligner, and then summarize the total depth and
coverage of each individual genome detected in the sample.

The workflow is as follows:

	1. Download any reference or sample data as needed
	2. Align raw reads against the reference database
	3. Keep all alignments for each read with amino acid similarity no more than 10% below than the best alignment
	4. Calculate the coverage and depth across each genome
	5. Copy the results in JSON format to the output location


This repository contains the information needed to build a Docker image, the preferred way
to run any bioinfomatic analysis. 

Only a single command is needed to run the analysis from start to finish:

```
map_viruses.py \
	--input <PATH_TO_INPUT_FILE> \
	--ref-db <REFERENCE_DATABASE> \
	--mapping <MAPPING_FILE> \
	--output_path <FILEPATH_FOR_OUTPUT> \
```

The `output_path` must end in `.json.gz`, as it will be formatted as a gzipped JSON file.

Read about additional parameters with `map_viruses.py --help`.


#### Input

Input files are FASTQ(.gz), and can be provided via local path, URL, or S3 key.


#### Reference database

Reference databases are indexed by DIAMOND, more details below.


#### Mapping

The mapping file links each reference protein to a genome of interest, more details below 
on how to build this file along as part of a reference database. 


#### Output

The output will be placed in a given directory in JSON.GZ format, including the following fields:

```
{
	"ref_db": STR,
	"input_path": STR,
	"logs": [STR, STR, ...],
	"time_elapsed": FLOAT,
	"input": STR,
	"total_reads": INT,
	"ref_db_url": STR,
	"output_folder": STR,
	"results": {
		"genomes": [
			{
		        "detected_proteins": 11,
		        "bitscore": 61.93033541343649,
		        "total_proteins": 11,
		        "pctid": 96.52610358441636,
		        "alen": 29.473899008492843,
		        "depth": 263.93339063171464,
		        "genome": "NC_001422.1",
		        "coverage": 0.9716373012462398,
		        "nreads": 20833,
		        "total_length": 2327
		    },
		    ...
		],
		"proteins": [
	        {
		        "definition": "internal scaffolding protein [Enterobacteria phage phiX174 sensu lato]",
		        "product": "internal scaffolding protein",
		        "bitscore": 62.25130237825595,
		        "alen": 29.160815402038505,
		        "taxonomy": "Viruses; ssDNA viruses; Microviridae; Bullavirinae; Phix174microvirus; unclassified Phix174microvirus",
		        "region": "5075..5386,NC_001422.1:1..51)",
		        "pctid": 96.95537938844848,
		        "locus_tag": "phiX174p03",
		        "taxid": 374840,
		        "length": 120,
		        "depth": 214.56666666666666,
		        "genome": "NC_001422.1",
		        "coverage": 1,
		        "nreads": 883,
		        "protein": "NP_040705",
		        "organism": "Enterobacteria phage phiX174 sensu lato"
	        },
	        ...
		]
	}
}
```

### Saving raw alignment files

If you would like to store the raw alignment files, use the `--keep-alignments` flag. This will copy
a file ending in ".sam.gz" to the output path that you specify, in addition to the other output files.


### Making a reference database

To make a reference database, simply create a FASTA file with the **protein** sequences for each virus,
and a tab-delimited text file linking each sequence to a specific genome. The only three columns that 
*must* be in that mapping file are `protein`, `length`, and `genome`.

Next, index the FASTA file with DIAMOND:

```
diamond makedb --in ref.fasta --db ref
```

Finally, when invoking `map_viruses.py`, use `--ref-db` to point to the DIAMOND indexed database, and 
`--mapping` to point to the tab-delimited text file linking the proteins to genomes.

**Automatically create a database**

You can automatically create a database by running the command `make_viral_db.py`. 
