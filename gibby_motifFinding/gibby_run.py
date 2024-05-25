import gibbs
from Bio import SeqIO
import numpy as np
import seqlogo
import os
import argparse

# Function to parse BED files
def parse_bed_file(peak_file, score_threshold):
    peaks = []
    with open(peak_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            score = float(fields[8])
            if score > score_threshold:
                peaks.append((chrom, start, end, score))
    return peaks

# Organizes the genome fasta file by chromosome so that we can later find sequences more quickly; MUCH FASTER
def get_genome_info(fasta_file):
    genome_info = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_info[record.id] = record.seq
    return genome_info

# Extract sequence from organized genome fasta file based on peak coordinates; FAST
def extract_sequence(chrom, start, end, genome_info):
    return str(genome_info[chrom][start:end])

# Get peak sequences
def get_peak_sequences(bed_file, genome_fasta_file, score_threshold):
    peaks = parse_bed_file(bed_file, score_threshold)  # Get peak coordinates
    genome_info = get_genome_info(genome_fasta_file)  # Organize genome fasta file for fast sequence extraction
    sequences = []  # Store peak sequences
    for peak in peaks:  # Extract sequences for each peak
        chrom, start, end, score = peak  # Unpack peak coordinates and score
        sequence = extract_sequence(chrom, start, end, genome_info)  # Extract sequence
        sequences.append(sequence.upper())  # Store sequence
    return sequences  # Return all sequences

nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

def GetPFM(sequences):
    pfm = np.zeros((4, len(sequences[0])))
    for sequence in sequences:
        for i in range(len(sequence)):
            if sequence[i] == "A":
                pfm[nucs['A'], i] += 1
            if sequence[i] == "T":
                pfm[nucs['T'], i] += 1
            if sequence[i] == "C":
                pfm[nucs['C'], i] += 1
            if sequence[i] == "G":
                pfm[nucs['G'], i] += 1
    return pfm

def find_background_freq(sequences):
    total = 0
    frequencies = {"A": 0, "C": 0, "G": 0, "T": 0}
    for sequence in sequences:
        for nucleotide in sequence:
            frequencies[nucleotide] += 1
            total += 1
    return {nucleotide: round(frequency / total, 2) for nucleotide, frequency in frequencies.items()}

def GetPWM(binding_sites, background_freqs=[0.25, 0.25, 0.25, 0.25]):
    """ Compute the PWM for a set of binding sites
    Parameters
    ----------
    binding_sites : list of str
        List of sequences 
    background_freqs: list of float
        Background frequency of A, C, G, T
    Returns
    -------
        pwm : 2d np.array
    Assumes all sequences have the same length
    """
    pwm = np.zeros((4, len(binding_sites[0])))
    pfm = GetPFM(binding_sites)
    pfm = pfm + 0.01 # Add pseudocount. Don't change this!
    # Compute pwm below
    # Note: np.sum(pfm[:,j]) will give the sum of counts for column j
    # Note: pfm[i,j]/np.sum(pfm[:,j]) gives p(i,j) (frequency of nucleotide i at position j)
    # your code here
    for i in range(len(pfm)):
        for j in range(len(pfm[i])):
                pwm[i,j] =  np.log2((pfm[i,j]/np.sum(pfm[:,j])) / background_freqs[i])
    return pwm

def main():
    parser = argparse.ArgumentParser(description="Process genomic data to extract peak sequences and compute PWMs.\n"
                                                "Usage: python gibby_run.py -b bed_file -g genome_fasta_file -s score_threshold -i iterations -k kmer_size")
    parser.add_argument('-b', '--bed_file', type=str, required=True, help="Input BED file with peak information.")
    parser.add_argument('-g', '--genome_fasta_file', type=str, required=True, help="Genome FASTA file.")
    parser.add_argument('-s', '--score_threshold', type=float, required=True, help="Score threshold for filtering peaks.")
    parser.add_argument('-k', '--kmer_size', type=int, required=True, help="Size of k-mers for motif finding.")
    parser.add_argument('-i', '--iterations', type=int, required=True, help="Number of iterations for Gibbs sampler.")

    args = parser.parse_args()

    dna = get_peak_sequences(args.bed_file, args.genome_fasta_file, args.score_threshold)

    parse_bed_file(args.bed_file, args.score_threshold)

    foo = gibbs.gibbs_sampler(dna, args.kmer_size, args.iterations)
    print("Gibbs Sampler Output: \n", foo)

    freq = find_background_freq(foo)
    print("Background Frequencies: \n", freq)

    pfm = GetPFM(foo)
    print("Position Frequency Matrix: \n", pfm)

    pwm = GetPWM(foo)
    print("Position Weight Matrix: \n", pwm)

    os.environ['PATH'] = '/opt/homebrew/bin:' + os.environ['PATH']

    seq_pwm = seqlogo.Pwm(pwm)
    seq_ppm = seqlogo.Ppm(seqlogo.pwm2ppm(seq_pwm))
    seqlogo.seqlogo(seq_ppm, ic_scale=False, format='png', size='large', filename='output_logo.png')

    print("Loaded and saved motif logo")

if __name__ == "__main__":
    main()


