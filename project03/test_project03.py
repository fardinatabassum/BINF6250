import random
import numpy as np
import bamnostic as bs
import seqlogo as sl
import shutil
sl.seqlogo.shutil = shutil

#import function for building sequence motif & idenfitying seqs matching to motif
from data_readers import *
from seq_ops import get_seq
from motif_ops import *

def get_kmers(seq, k):
    kmer_seqs = []
    seq_rev = reverse_complement(seq)
    for ind in range(0, len(seq) - (k-1)):
        kmer_seqs.append(seq[ind:ind+k])
        kmer_seqs.append(seq_rev[ind:ind+k])

    return kmer_seqs

def log2_weight_transform(scores):
    for ind, score in enumerate(scores):
        scores[ind] = 2**score

    return scores


def GibbsMotifFinder (seqs, k, seed=None):
    '''
    Function to find a pfm from a list of strings using a Gibbs sampler
    
    Args: 
        seqs (str list): a list of sequences, not necessarily in same lengths
        k (int): the length of motif to find
        seed (int, default=None): seed for np.random

    Returns:
        pfm (numpy array): dimensions are 4xlength
    '''
    monocharacter_seqs = {"A"*k, "C"*k, "T"*k, "G"*k}
    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    # Get initial motifs for pfm
    motifs = []
    for seq in seqs:
        randint = random.randint(0, len(seq) - k)
        motifs.append(seq[randint:randint+k])

    for _ in range(50000):
        rand_seq_ind = random.randint(0, len(seqs)-1)
        rand_seq = seqs[rand_seq_ind]

        kmer_seqs = get_kmers(rand_seq, k)

        pfm = build_pfm(motifs, k)
        pwm = build_pwm(pfm)

        valid_kmers = []
        kmer_scores = []

        for kmer in kmer_seqs:
            if kmer in monocharacter_seqs:
                continue
            kmer_scores.append(score_kmer(kmer, pwm))
            valid_kmers.append(kmer)

        #print(kmer_scores)
        #print("W"*80)
        weights = log2_weight_transform(kmer_scores)
        #print(weights)
        

        new_motif = random.choices(valid_kmers, weights=weights)[0]
        motifs[rand_seq_ind] = new_motif
        if _ % 2000 == 0:
            print(len(kmer_seqs))
            print(f"Iteration {_} IC Score: {pfm_ic(pfm)}", end="\n"*2)
    print(new_motif)
    return pfm


if __name__ == "__main__":
    """bam_path = "data/SRR9090854.subsampled_5pct.bam"
    print("Reading in sequences")
    seqs = [read.seq for read in bs.AlignmentFile(bam_path)]
    with open("data/ChIP_Subset.txt", mode = 'w', encoding='utf-8') as txtout:
        for seq in seqs[:10000]:
            print(seq, file=txtout)"""
    #with open("data/ChIP_Subset.txt", mode='r', encoding='utf-8') as data_file:
    #    seqs = [read.rstrip() for read in data_file]
    with open("data/peaks100.fa", mode='r', encoding='utf-8') as data_file:
        seqs = [read.rstrip().upper() for read in data_file if not read[0] == ">" and not "N" in read]
    print(len(seqs))
    """print(seqs[:5])
    for ind, n in enumerate(seqs):
        if "N" in n:
            print("Oh Boy")
            print(n)
    
    print(seqs[107])"""
            

    # Run the gibbs sampler:
    promoter_pfm = GibbsMotifFinder(seqs, 10)
    print(promoter_pfm)

    # Plot the final pfm that is generated: 
    #sl.seqlogo(sl.CompletePm(pfm = promoter_pfm.T), format='pdf')

