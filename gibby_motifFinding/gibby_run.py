import argparse
from Bio import SeqIO
import numpy as np
import seqlogo
import os
import random
from random import choices

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
    pwm = np.zeros((4, len(binding_sites[0])))
    pfm = GetPFM(binding_sites)
    pfm = pfm + 0.01  # Add pseudocount. Don’t change this!
    for i in range(len(pfm)):
        for j in range(len(pfm[i])):
            pwm[i, j] = np.log2((pfm[i, j] / np.sum(pfm[:, j])) / background_freqs[i])
    return pwm

def gibbs_sampler(dna, k, n):
    """Implements the GibbsSampling algorithm for motif finding."""
    best_motifs = []
    best_score = float('inf')
    t = len(dna)
    all_motifs_scores = []

    for _ in range(1000):
        motifs = randomKmers(dna, k)  # random Kmer Patterns and their start positions
        current_best_motifs = motifs[:]  # copy motifs for initialization
        current_best_score = Score(current_best_motifs, k)  # score the motifs and initialize the best score
        for _ in range(n):
            ignoredKmerIndex = chooseRandomKmerIndex(dna)  # choose random i
            profile = Profile(current_best_motifs[:ignoredKmerIndex] + current_best_motifs[ignoredKmerIndex + 1:], k)  # make a profile from the kmers excluding the ignored Kmer
            motifProbability = Motifs(profile, dna[ignoredKmerIndex], k)  # find motif probability from the ignoredDNA
            probableKmers = list(motifProbability.keys())  # kmers in ignoredDNA as they are the keys in the profile
            probabilities = list(motifProbability.values())  # probabilities of those kmers using profile
            chosenKmer = choices(probableKmers, probabilities)  # choose a kmer using the probabilities
            chosenKmer = "".join(chosenKmer)  # join because it comes out in list format
            current_best_motifs[ignoredKmerIndex] = chosenKmer  # swap out the ignored motif with the chosenKmer
            score = Score(current_best_motifs, k)
            all_motifs_scores.append((chosenKmer, score))
            if score < current_best_score:  # score comparison
                current_best_score = score  # If the score is better replace current_best_score
            if score < best_score:
                best_motifs = current_best_motifs[:]  # replace best as well for both. Leave this because it is the final answer
                best_score = score

    return best_motifs

def randomKmers(dna, k):
    randomK = []
    for strings in dna:
        if len(strings) >= k:  # Ensure the string is at least as long as k
            start_pos = random.randint(0, len(strings) - k)  # randomly choose kmer start position
            pattern = strings[start_pos:start_pos + k]  # the kmer
            randomK.append(pattern)
    return randomK

def chooseRandomKmerIndex(dna):
    t = len(dna)
    index = random.randint(0, t - 1)  # choose a random index
    return index

def Profile(randomKmers, k):
    profile = [{'A': 1, 'C': 1, 'G': 1, 'T': 1} for _ in range(k)]  # make dictionaries in list for every k START AT 1 FOR Pseudocounts
    for i in range(k):
        for string in randomKmers:
            base = string[i]
            profile[i][base] += 1  # add one to the nucleotide count based on the letter
    for dicts in profile:
        total = sum(dicts.values())  # sum so that divide by the total to get probability
        for nucleotide in dicts:
            dicts[nucleotide] /= total  # probability

    return profile

def Motifs(profile, string, k):
    motifProbability = {}
    for i in range(len(string) - k + 1):
        kmer = string[i:i + k]  # set window
        probability = 1
        for j in range(k):
            base = kmer[j]  # base at j
            probability *= profile[j][base]  # probability at position j given the base
        motifProbability[kmer] = probability  # put kmer and probability in dict so that we can split later for randomly "choosing" given probability
    return motifProbability

def Score(mostProbable, k):
    score = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for _ in range(k)]
    for i in range(k):
        for string in mostProbable:  # for every kmer in mostProbable
            base = string[i]  # take the base
            score[i][base] += 1  # this makes the count for each nucleotide in the i-th position
    sumList = []
    for i in range(k):
        dicts = score[i]  # all the dictionaries in the score
        mostFreq = max(dicts.values())  # most frequent base that is seen in that position
        total = len(mostProbable)  # total is the number of how many motifs there are
        indexScore = total - mostFreq  # total - most frequent base leads to the indexScore
        sumList.append(indexScore)
    return sum(sumList)  # sum it all up

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

    foo = gibbs_sampler(dna, args.kmer_size, args.iterations)
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



