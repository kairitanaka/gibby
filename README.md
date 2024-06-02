![Python Version](https://img.shields.io/badge/python-%3E%3D3.6-red)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
![Open Issues](https://img.shields.io/github/issues/kairitanaka/gibby?color=yellow)
![Contributors](https://img.shields.io/github/contributors/kairitanaka/gibby?color=green)

# Gibby (v0.1.0)

`Gibby` is a Python package designed to find the motif of a transcription factor de novo based on ChIP-seq data. Given a genome and peak file, the tool extracts sequences from the genome based on the peak file, and runs Gibbs Sampling on those peak sequences. During the sampling process, kmers across the peak sequences are compared and scored. The output of Gibbs sampling is a list potential motifs of length k; the position frequency and position weight matrices of these potential motifs are created and saved as a text file. The position weight matrix is visualized using seqLogo to observe which motif was most strongly conserved among the peak regions of the transcription factor. 

- [What's the input?](#input)
- [What's the output?](#output)
- [What is Gibbs Sampling?](#what-is-gibbs-sampling)
- [All Features](#features)
- [How do I install Gibby?](#installation)
- [How do I use Gibby?](#usage)
- [Examples](#examples)

## Input 
`Gibby` can take HOMER `peaks.txt` files and `.bed` files generated from ChIP-seq data as input. In addition, we require the appropriate genome assembly (which was used to generate the peak file) in `FASTA` format. 

## Output
Three files are generated: `PFM.txt`, `PWM.txt`, and `motif.png`. These files include the position frequency matrix and position weight matrix generated by gibbs sampling, and a visualization of the sequence motif showing which motif was most strongly conserved among the peak regions. 

<p align="center">
<img src=https://github.com/kairitanaka/CSE_185_finalProject/assets/86521451/f60bc209-1975-4fe0-9844-40f2f0de3f25)>
</p>

# What is Gibbs Sampling?
Gibbs Sampling is a statistical method used to estimate the distribution of variables when direct sampling is difficult. It is particularly useful in motif finding where we want to identify common patterns (motifs) in a set of sequences. 

## For example ...
Imagine you are trying to perfect a secret recipe, but you don't have all the ingredients at once. You start by randomly choosing some ingredients and proportions, then you taste the result. Based on how good or bad it tastes, you keep some ingredients, change others, and try again. Each iteration helps you understand what works and what doesn't, gradually leading you to the best recipe.

In the context of genomic sequences, Gibbs Sampling helps us identify common patterns (motifs) by iteratively refining our guesses based on the sequences we have. Each iteration, which involves some randomness, helps to gradually reveal the underlying motifs more accurately.

## Process for Gibbs Sampling

**The problem we are trying to solve here is:**\
Given `S` sequences, find the most mutually similar subsequences of length `k` from each sequence

In order to tackle this problem it is crucial to look at the entire statistical landscape by sampling every single sequence and seeing if we can converge to a minima that is the optimal or somewhere extremely close to the optimum. 

Randomly choose a starting position for the subsequence of length `k` in each of the `S` sequences.
For each iteration, leave out one sequence, say sequence `s'`.

<p align="center">
  <img src=https://github.com/kairitanaka/gibby/assets/64274901/b48c225e-da4e-4d6a-b274-aaf5bd5ce36f width="50%" alt="Image 1 Description">
</p>

Using the remaining `S-1` sequences, build a position-specific scoring matrix (PSSM) or profile matrix. This matrix represents the frequency of each nucleotide at each position of the subsequence. USE PSEUDOCOUNTS!!!

<p align="center">
  <img src=https://github.com/kairitanaka/gibby/assets/64274901/87bd1dbc-7069-4840-9144-a51a5f151936 width="50%" alt="Image 1 Description">
</p>


Calculate the probability of every possible subsequence of length `k` in the left-out sequence `s'` using the profile matrix. This involves calculating the likelihood of the subsequence given the profile and normalizing it to get a probability distribution.

<p align="center">
  <img src=https://github.com/kairitanaka/gibby/assets/64274901/b9a044a8-b4fb-489e-bc7c-2c7bf0359d9e width="50%" alt="Image 1 Description">
</p>

$$Random(\frac{2}{8^4}, \frac{2}{8^4}, \frac{72}{8^4}, \frac{24}{8^4}, \frac{8}{8^4},\frac{4}{8^4}, \frac{1}{8^4}) = Random(\frac{2}{113},\frac{2}{113},\frac{72}{113},\frac{24}{113},\frac{8}{113},\frac{4}{113},\frac{1}{113})$$

Sample a new position for the subsequence in sequence `s'` according to the probability distribution obtained in the previous step. This new position replaces the old position for sequence `s'`. In this case we chose ACCT. 

<p align="center">
  <img src=https://github.com/kairitanaka/gibby/assets/64274901/89a325c2-edde-4aba-b040-c2224206c672 width="50%" alt="Image 1 Description">
</p>


After this we will score the motifs. If the score is better than previous iteration we will keep those set of motifs so that we are always progressing in the correct direction. Then we will repeat for multiple iterations!

We have seen that in around 500 - 1000 iterations the positions of the subsequences have stabilized across iterations. However, this may take some testing over 2~5 runs based on your dataset. 

Reference: [bioinformatics algorithms an active learning approach](https://www.bioinformaticsalgorithms.org/) and BIMM 181

## TL;DR: A short explanation of Gibbs Sampling
Steps for Gibbs Sampling:

Another resource you can use is this [video](https://www.youtube.com/watch?v=MP6O_Z2AUDU) from our beloved professor from UCSD,  Dr.Pavel. 

# Features

- Utilize Gibbs sampling to **identify motifs** within the peak regions for a given transcription factor.
- Compute Position Frequency/Weight Matrices (PFMs, PWMs)

# Installation

You can install the package via pip:

```
pip install git+https://github.com/kairitanaka/gibby.git
```
and verify by running:
```
gibby -h
```
If you get an error that the command is not found, make sure the directory ~/.local/bin is included in your $PATH environment variable. Or consider adding the directory to your `$PATH` by running:
```
export PATH=$PATH:$HOME/.local/bin
```
This will allow you to simply type `gibby` to run the tool. You will have to repeat this for every new terminal session.

If you come across the error:
```
error: cannot find command 'git'
```
Please make sure to install Git. You can find installation instructions [here](https://github.com/git-guides/install-git).


# Basic Usage

`gibby`, from Gibby (ver 0.1.0), utilizes Gibbs Sampling to find potential motifs that are in peak regions of the genome. The potential motifs are used to generate a position frequency matrix (PFM), a position weight matrix (PWM), and a motif logo based on these matrices. The options appear as below:

```
gibby [-h] -p PEAK_FILE -t PEAK_FILE_TYPE -g GENOME_FASTA_FILE [-s SCORE_THRESHOLD] [-k KMER_SIZE] [-i ITERATIONS]
```
### Required arguments
Assuming you have successfully installed Gibby, running the tool is a fairly simple task. First, it's a good idea to run `gibby -h` to see what options are available. You will notice that Gibby will always require three arguments to be passed: the `PEAK_FILE`, the `PEAK_FILE_TYPE` ("bed" or "homer"), and the `GENOME_FASTA_FILE`. In addition, there are several optional arguments such as `SCORE_THRESHOLD`, `KMER_SIZE`, and `ITERATIONS`. 

### Optional arguments
`SCORE_THRESHOLD` represents the minimum score that the tool will use to filter out low quality peaks. Generally, if you know the general length of the motif that you are looking for, you can specify `KMER_SIZE` to the general length + 5 (to add some leeway). `ITERATIONS` is the number of times you want to run Gibbs Sampling. Running more iterations results in better detection of motifs in exchange for the increasing time it takes to finish the task. We have set default values for these variables which strike a good balance between speed and performance.

## Command Line Arguments

- `-p`, `--peak_file`: Input BED/Homer file with peak information. (required)
- `-t`, `--peak_file_type`: Specify peak file type: 'bed' or 'homer' file. (required)
- `-g`, `--genome_fasta_file`: Genome FASTA file. (required)
- `-s`, `--score_threshold`: Score threshold for filtering peaks. (optional)
- `-k`, `--kmer_size`: Size of k-mers for motif finding. (optional)
- `-i`, `--iterations`: Number of iterations for Gibbs sampler. (optional)


# Examples

## OCT4 (Pou5f1) 
In this example, we will be using the same files that were used in Lab 5 of the CSE 185 course. The files you will need are the following: (1) `peaks.txt` which is the peak file generated by HOMER for OCT4 (which we have also included in our repo) and (2) the `GRCm38.fa` Mus musculus genome assembly in FASTA format which should be available in the CSE 185 public genomes folder. 

In this case, we have a HOMER peak file. In addition, suppose we want to choose 75 as the score threshold to filter out low-quality peaks. Running the tool would look like this: 

```
gibby -p ~/lab5/tagdirs/Oct4/peaks.txt -t homer -g ~/public/genomes/GRCm38.fa -s 75
```

You will want to make sure that you have the correct paths for the peak file and the genome file. Running the command, the tool will take some time to fully run. After finishing, you will want to take a look at the generated motif logo which visualizes which motif was most conserved among the peak regions. In our case, we got: 
<p align="center">
<img src=https://github.com/kairitanaka/CSE_185_finalProject/assets/86521451/437e5a4a-e57f-4b51-be1b-27c271150b64)>
</p>

and sometimes its reverse complement:
<p align="center">
<img src=https://github.com/kairitanaka/CSE_185_finalProject/assets/86521451/bd12cc26-be4a-4e9e-aa6b-d68610c6b8e1)>
</p>

Since Gibbs Sampling is a stochastic process, your graphs will look different from our results. However, what *should* be similar are the large letters (nucleotides) stacked at the top and their relative positioning to one another. These large sets of letters represent the motif that Gibby discovered; all the other smaller nucleotides represent "noise" or nucleotides that were not as strongly conserved among peak regions for OCT4. You may get the motif or its reverse complement. 

Compare your results with the published motifs for OCT4: 
<p align="center">
<img src=https://github.com/kairitanaka/gibby/assets/86521451/af7635bd-5fdc-4d41-8594-a5e0b9108157)>
</p>

Reverse complement:
<p align="center">
<img src=https://github.com/kairitanaka/gibby/assets/86521451/4daf935c-9153-45e5-9bf7-a1814e91964a)>
</p>

Feel free to test the tool with the other two transcription factors (KLF4 and SOX2) used in Lab 5. Since they use the same genome assembly, you just need to change the HOMER peak file!

## ZNF24
In this example, we will be using a ChIP-seq dataset that is for the transcription factor ZNF24. The bed file we used is from ENCODE: `https://www.encodeproject.org/files/ENCFF664TYB/`. The GRCh38 genome assembly can be downloaded from the UCSC Genome Browser: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/`. 

Since we are using a `.bed` peak file, we can configure the argument accordingly and use default settings for the other arguments:
```
gibby -p ENCFF664TYB.bed -t bed -g hg38.fa
```
Again, make sure the paths to the files are correct. Below we share the motif result we got, and the published motif for ZNF24. Gibby can return the motif and sometimes its reverse complement due to its stochasticity. 

#### Gibby Results:
Forward:
<p align="center">
<img src=https://github.com/kairitanaka/gibby/assets/86521451/3628e112-6fe1-4e91-8901-4de2bae589cf)>
</p>

Reverse Complement:
<p align="center">
<img src=https://github.com/kairitanaka/gibby/assets/86521451/352455f1-9ef7-462e-a0a8-a016c4057ea8)>
</p>

#### Published motif:
Forward:
<p align="center">
<img src=https://github.com/kairitanaka/gibby/assets/86521451/ea6efb5b-c0f8-4c88-b1b5-67394f5e8f5a)>
</p>

Reverse Complement:
<p align="center">
<img src=https://github.com/kairitanaka/gibby/assets/86521451/9c680cd2-17a5-47d4-a18b-44628ab7e04b>
</p>

## That's it! Hopefully they look similar!
## Authors 
Joe Hwang (j8hwang@ucsd.edu) & Kairi Tanaka (ktanaka@ucsd.edu)

