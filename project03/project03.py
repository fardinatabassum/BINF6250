import random
import numpy as np
import bamnostic as bs
import seqlogo as sl

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

        weights = log2_weight_transform(kmer_scores)
        

        new_motif = random.choices(valid_kmers, weights=weights)[0]
        motifs[rand_seq_ind] = new_motif
        if _ % 2000 == 0:
            print(len(kmer_seqs))
            print(f"Iteration {_} IC Score: {pfm_ic(pfm)}", end="\n"*2)
    print(new_motif)
    return pfm


if __name__ == "__main__":
    seq_file="data/GCF_000009045.1_ASM904v1_genomic.fna"
    gff_file="data/GCF_000009045.1_ASM904v1_genomic.gff"

    seqs = []
    for name, seq in get_fasta(seq_file):
        for gff_entry in get_gff(gff_file):
            if gff_entry.type == 'CDS':
                promoter_seq = get_seq(seq, gff_entry.start, gff_entry.end, gff_entry.strand, 50)
                if "AGGAGG" in promoter_seq:
                    seqs.append(promoter_seq)

    # Run the gibbs sampler:
    promoter_pfm = GibbsMotifFinder(seqs, 10)
    print(promoter_pfm)

    # Plot the final pfm that is generated: 
    #sl.seqlogo(sl.CompletePm(pfm = promoter_pfm.T), format='pdf')
