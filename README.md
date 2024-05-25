## What is Gibbs Sampling?
Gibbs Sampling is a statistical method used to estimate the distribution of variables when direct sampling is difficult. It is particularly useful in motif finding where we want to identify common patterns (motifs) in a set of sequences. 
![gibbs](https://github.com/kairitanaka/CSE_185_finalProject/assets/64274901/4d090d70-eebb-4224-9396-825bc4142684)

# For example ...
Imagine you are trying to perfect a secret recipe, but you don't have all the ingredients at once. You start by randomly choosing some ingredients and proportions, then you taste the result. Based on how good or bad it tastes, you keep some ingredients, change others, and try again. Each iteration helps you understand what works and what doesn't, gradually leading you to the best recipe.

In the context of genomic sequences, Gibbs Sampling helps us identify common patterns (motifs) by iteratively refining our guesses based on the sequences we have. Each iteration, which involves some randomness, helps to gradually reveal the underlying motifs more accurately.

More images and text to come explaining gibbs ... 


# gibby_motifFinding

`gibby_run` is a Python package designed to process genomic data, extract peak sequences, compute Position Weight Matrices (PWMs), and compare scores in order to converge to the best pattern. This tool is particularly useful for motif finding and other genomic analyses.

## Input 
It takes a HOMER peaks.txt or a bed peak file. The output will be a motif logo in the directory you run it in!

## Features

- **Extract Peak Sequences**: Retrieve sequences from genomic data based on peak information.
- **Compute Position Weight Matrices (PWMs)**: Generate PWMs from the extracted sequences.
- **Gibbs Sampling for Motif Finding**: Utilize Gibbs sampling to identify motifs within the sequences.

## Installation

You can install the package via pip:

```
pip install git+https://github.com/kairitanaka/CSE_185_finalProject.git
```
Once you run it it will be helpful to run:
```
gibby_motifFinding -h
```
for all the options


## Usage

The `gibby_motifFinding` script processes genomic data to extract peak sequences and compute PWMs using Gibbs sampling. The script can be run with the following command:

```
gibby_motifFinding -b bed_file -g genome_fasta_file -s score_threshold -i iterations -k kmer_size
```


### Command Line Arguments

- `-b`, `--bed_file`: Input BED file with peak information. (required)
- `-g`, `--genome_fasta_file`: Genome FASTA file. (required)
- `-s`, `--score_threshold`: Score threshold for filtering peaks. (required)
- `-k`, `--kmer_size`: Size of k-mers for motif finding. (required)
- `-i`, `--iterations`: Number of iterations for Gibbs sampler. (required)



