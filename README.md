# Map Viruses
Align short reads against viral genomes in protein space


Researchers frequently want to detect viral genomes from short read metagenomic datasets. 
My preferred way of doing this is to align reads against those viral genomes in amino acid
space (akin to BLASTx) using the DIAMOND aligner, and then summarize the total depth and
coverage of each individual genome detected in the sample.

The workflow is as follows:

	1. Download any reference or sample data as needed
	2. Align raw reads against the reference database
	3. Keep all alignments for each read with amino acid similarity no more than 10% below than the best alignment
	4. Calculate the coverage and depth across each genome
	5. Copy the results in JSON format to the output directory


This repository contains the information needed to build a Docker image, the preferred way
to run any bioinfomatic analysis. Only a single command is needed to run the analysis from
start to finish:

```
map_viruses.py \
	--input <PATH_TO_INPUT_FILE> \
	--ref-db <REFERENCE_DATABASE> \
	--mapping <MAPPING_FILE> \
	--output <DIRECTORY_FOR_OUTPUT> \
```

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
	"input": STR,
	"results": [
		{
			"genome": STR,
			"species": STR,
			"coverage": FLOAT,
			"depth": FLOAT,
			"nreads": INT,
		}
	]
}
```


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
