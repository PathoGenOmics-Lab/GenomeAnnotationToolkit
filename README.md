# GenomeAnnotationToolkit
This repository contains useful scripts for parsing genomic annotations.

## Scripts available
1. [gbk2tsv.py: GenBank to TSV Converter](#gbk2tsvpy-genbank-to-tsv-converter)

# gbk2tsv.py: GenBank to TSV Converter

This Python script is designed to convert GenBank files into tab-separated value (TSV) files. It is particularly useful for bioinformaticians and researchers studying gene sequences. 
Modified from [this repository](https://github.com/wanyuac/BINF_toolkit/blob/master/gbk2tsv.py).

## Features

- Conversion of GenBank files into TSV format
- Selection of specific features to include (CDS, tRNA, rRNA etc.)
- Option to include nucleotide sequences of features
- Option to include protein sequences of Coding Sequences (CDS)

## Requirements

- Python 3.x
- Biopython
- argparse

## Usage
```
python3 gbk2tsv.py --gbk file.gbk --outdir ./output_directory --features CDS,tRNA,rRNA --nucleotides --protein
```

You can also process multiple GenBank files at once:
python3 gbk2tsv.py --gbk $(ls *.gbk) --outdir . --features "CDS,rRNA,tRNA"
Or specify individual GenBank files:

```
python3 gbk2tsv.py --gbk file1.gbk file2.gbk --outdir . --features CDS --nucleotides
```

## Arguments

- `-g`, `--gbk`: Specify GenBank files as input
- `-o`, `--outdir`: Define the directory for output files
- `-f`, `--features`: Features to include, separated by commas
- `-n`, `--nucleotides`: Enable this to include nucleotide sequences of features
- `-p`, `--protein`: Enable this to include protein sequences of CDS

## Logging

The script generates a log file with detailed information about the processing of each file. The log file is named `gbk2tsv_YYYY-MM-DD_HH-MM-SS.log`, with the date and time reflecting when the script was run.
