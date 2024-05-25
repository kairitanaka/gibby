# gibby_motifFinding

`gibby_run` is a Python package designed to process genomic data, extract peak sequences, and compute Position Weight Matrices (PWMs). This tool is particularly useful for motif finding and other genomic analyses.

## Features

- **Extract Peak Sequences**: Retrieve sequences from genomic data based on peak information.
- **Compute Position Weight Matrices (PWMs)**: Generate PWMs from the extracted sequences.
- **Gibbs Sampling for Motif Finding**: Utilize Gibbs sampling to identify motifs within the sequences.

## Installation

You can install the package via pip:

```
python gibby_run.py -b bed_file -g genome_fasta_file -s score_threshold -i iterations -k kmer_size
```


## Usage

The `gibby_run` script processes genomic data to extract peak sequences and compute PWMs using Gibbs sampling. The script can be run with the following command:

```
python gibby_run.py -b bed_file -g genome_fasta_file -s score_threshold -i iterations -k kmer_size
```


### Command Line Arguments

- `-b`, `--bed_file`: Input BED file with peak information. (required)
- `-g`, `--genome_fasta_file`: Genome FASTA file. (required)
- `-s`, `--score_threshold`: Score threshold for filtering peaks. (required)
- `-k`, `--kmer_size`: Size of k-mers for motif finding. (required)
- `-i`, `--iterations`: Number of iterations for Gibbs sampler. (required)



