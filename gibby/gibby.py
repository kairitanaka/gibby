import argparse
from Bio import SeqIO
import numpy as np
import seqlogo
import os
import random
from random import choices
from tqdm.auto import tqdm
import threading
import time
import sys
from termcolor import colored
from pyfiglet import Figlet

# Custom error class to handle exceptions
class CustomError(Exception):
    pass

# Function to parse BED files
def parse_bed_file(peak_file, file_type, score_threshold):
    try:
        peaks = []
        # Parse BED file
        if file_type == "bed":
            with open(peak_file, 'r') as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    score = float(fields[4])
                    if score > score_threshold:
                        peaks.append((chrom, start, end, score))
            return peaks
        # Parse HOMER file
        elif file_type == "homer":
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
        # Raise error if file type is not 'bed' or 'homer'
        else:
            raise CustomError("Invalid file type. Please specify 'bed' or 'homer'.")
    except Exception as e:
        raise CustomError(f"\n ERROR PARSING FILE: {str(e)}")
    return peaks

# Organizes the genome fasta file by chromosome so that we can later find sequences more quickly; MUCH FASTER
def get_genome_info(fasta_file):
    try:
        genome_info = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            genome_info[record.id] = record.seq
        return genome_info
    except Exception as e:
        raise CustomError(f"\n ERROR READING GENOME FASTA FILE: {str(e)}")

# Extract sequence from organized genome fasta file based on peak coordinates.
def extract_sequence(chrom, start, end, genome_info):
    try:
        return str(genome_info[chrom][start:end])
    except Exception as e:
        raise CustomError(f"\n ERROR EXTRACTING SEQUENCE: {str(e)}")

# Get peak sequences
def get_peak_sequences(bed_file, file_type, genome_fasta_file, score_threshold=None):
    try:
        if file_type == "bed" and score_threshold is None:
            score_threshold = 999 # Default score threshold for BED files
        elif file_type == "homer" and score_threshold is None:
            score_threshold = 60 # Default score threshold for HOMER files
        peaks = parse_bed_file(bed_file, file_type, score_threshold)  # Get peak coordinates
        genome_info = get_genome_info(genome_fasta_file)  # Organize genome fasta file for fast sequence extraction
        sequences = []  # Store peak sequences
        for peak in peaks:  # Extract sequences for each peak
            chrom, start, end, score = peak  # Unpack peak coordinates and score
            sequence = extract_sequence(chrom, start, end, genome_info)  # Extract sequence
            sequences.append(sequence.upper())  # Store sequence
        return sequences  # Return all sequences
    except Exception as e: # Catch any exceptions
        print(f"{e}")

nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

# Get Position Frequency Matrix
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

# Get Position Weight Matrix
def GetPWM(binding_sites):
    pwm = np.zeros((4, len(binding_sites[0])))
    pfm = GetPFM(binding_sites)
    pfm = pfm + 0.01  # Add pseudocount. Don’t change this!
    for i in range(len(pfm)):
        for j in range(len(pfm[i])):
            pwm[i, j] = np.log2((pfm[i, j] / np.sum(pfm[:, j])) / 0.25)
    return pwm

# Gibbs Sampling Algorithm
def gibbs_sampler(dna, k, n):
    """Implements the GibbsSampling algorithm for motif finding."""
    if k is None: 
        k = 20 # default kmer size
    if n is None:
        n = 500 # default number of iterations
    best_motifs = []
    best_score = float('inf')
    t = len(dna)
    all_motifs_scores = []

    for i in tqdm(range(1000)):
        motifs = randomKmers(dna, k)  # random Kmer Patterns and their start positions
        current_best_motifs = motifs[:]  # copy motifs for initialization
        current_best_score = Score(current_best_motifs, k)  # score the motifs and initialize the best score
        for j in range(n):
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

def log_message(message):
    print(message)

def main():
    ascii_art = r"""
   _____   _   _       _             
  / ____| (_) | |     | |            
 | |  __   _  | |__   | |__    _   _ 
 | | |_ | | | | '_ \  | '_ \  | | | |
 | |__| | | | | |_) | | |_) | | |_| |
  \_____| |_| |_.__/  |_.__/   \__, |
                                __/ |
                               |___/            
    """
    explanation = ("\n Gibby (ver 0.1.0) utilizes Gibbs Sampling to find potential motifs that are in peak regions of the genome."
                   "The potential motifs are used to generate a position frequency matrix (PFM), a position weight matrix (PWM),"
                   "and a motif logo based on these matrices. \n")
    
    parser = argparse.ArgumentParser(description=ascii_art + explanation, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--peak_file', type=str, required=True, help="Peak file in BED or HOMER format.")
    parser.add_argument('-t', '--peak_file_type', type=str, required=True, help="File type: 'bed' or 'homer' file.")
    parser.add_argument('-g', '--genome_fasta_file', type=str, required=True, help="Genome FASTA file.")
    parser.add_argument('-s', '--score_threshold', type=float, required=False, default=None, help="Score threshold for filtering peaks.")
    parser.add_argument('-k', '--kmer_size', type=int, required=False, default=20, help="Size of k-mers for motif finding.")
    parser.add_argument('-i', '--iterations', type=int, required=False, default=500, help="Number of iterations for Gibbs sampler.")
    args = parser.parse_args()
    
    # Add homebrew path to PATH
    os.environ['PATH'] = '/opt/homebrew/bin:' + os.environ['PATH']

    # Function to write matrix to file
    def write_matrix_to_file(array, filename):
        np.savetxt(filename, array, fmt='%s')

    # Function to display loading dots animation
    def loading_dots():
        while not stop_loading:
            for i in range(4):
                if stop_loading:
                    break
                print(f'\rLoading peak sequences{"." * i}{" " * (3 - i)}', end='')
                time.sleep(0.5)

    # Main code
    try:
        stop_loading = False 
        loading_thread = threading.Thread(target=loading_dots)
        loading_thread.start()
        # Get peak sequences
        dna = get_peak_sequences(args.peak_file, args.peak_file_type, args.genome_fasta_file, args.score_threshold)
        stop_loading = True
        loading_thread.join()
    except CustomError as e:
        print(f"Error: {str(e)}")
        stop_loading = True
        loading_thread.join()
        sys.exit(1)
    
    if dna is None:
        print("Failure to find peak sequences. Exiting...")
        sys.exit(1)

    print("\n")
    print(len(dna), "peak sequences loaded. \n")
    print("Now running gibbs sampler...")

    # Run Gibbs Sampler
    foo = gibbs_sampler(dna[:100], args.kmer_size, args.iterations)

    print("Finished! \n")
    print("Writing potential binding motif to potential_motifs.txt...")

    # Write potential motifs to file
    write_matrix_to_file(foo, 'gibbs_potential_bind.txt')

    print("Done! \n")

    # Generate PFM, PWM, and motif logo
    pfm = GetPFM(foo)
    print("Writing position frequency matrix to PFM.txt...")
    write_matrix_to_file(pfm, 'PFM.txt')
    log_message(colored("Done writing PFM. \n", 'green'))
    pwm = GetPWM(foo)
    print("Writing position weight matrix to PWM.txt...")
    write_matrix_to_file(pwm, 'PWM.txt')
    log_message(colored("Done writing PWM. \n", 'green'))
    print("Generating motif logo...")
    seq_pwm = seqlogo.Pwm(pwm)
    seq_ppm = seqlogo.Ppm(seqlogo.pwm2ppm(seq_pwm))
    seqlogo.seqlogo(seq_ppm, ic_scale=False, format='png', size='large', filename='motif.png')

    log_message(colored("Successfully loaded and saved motif logo to motif.png \n", 'green'))

    completion_banner = Figlet(font='slant').renderText('Task Completed!')
    print(colored(completion_banner, 'cyan'))

if __name__ == "__main__":
    main()
