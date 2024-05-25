import sys
import random
from random import choices
from collections import Counter

def gibbs_sampler(dna: list[str], k: int, n: int) -> list[tuple[str, int]]:
    """Implements the GibbsSampling algorithm for motif finding."""
    best_motifs = []
    best_score = float('inf')
    t = len(dna)
    all_motifs_scores = []

    # We no longer calculate background frequencies since we're not using them

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

def randomKmers(dna: list[str], k: int) -> list:
    randomK = []
    for strings in dna:
        if len(strings) >= k:  # Ensure the string is at least as long as k
            start_pos = random.randint(0, len(strings) - k)  # randomly choose kmer start position
            pattern = strings[start_pos:start_pos + k]  # the kmer
            randomK.append(pattern)
    return randomK

def chooseRandomKmerIndex(dna: list) -> int:
    t = len(dna)
    index = random.randint(0, t - 1)  # choose a random index
    return index

def Profile(randomKmers: list, k: int):
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

def Motifs(profile, string, k) -> dict:
    motifProbability = {}
    for i in range(len(string) - k + 1):
        kmer = string[i:i + k]  # set window
        probability = 1
        for j in range(k):
            base = kmer[j]  # base at j
            probability *= profile[j][base]  # probability at position j given the base
        motifProbability[kmer] = probability  # put kmer and probability in dict so that we can split later for randomly "choosing" given probability
    return motifProbability

def Score(mostProbable: list, k: int):
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



