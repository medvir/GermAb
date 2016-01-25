#!/usr/local/bin/python3.4
'''
Trims random nucleotides used at the beginning of read 1 in Illumina amplicon sequencing.
Trims forward and reverse primers.
Allows for a certain number of mismatches between primer and nucleotide sequences (not conting wobble bases).
Arguments are:
1. sequences (fasta)
2. forward primers (fasta)
3. reverse primers (fasta)
Output prints to standard output.
'''

import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
 

dna_code = {
     'A': set(['A']),
     'C': set(['C']),
     'G': set(['G']),
     'T': set(['T']),
     
     'R': set(['G', 'A']),
     'Y': set(['T', 'C']),
     'M': set(['A', 'C']),
     'K': set(['G', 'T']),
     'S': set(['G', 'C']),
     'W': set(['A', 'T']),
     
     'H': set(['A', 'C', 'T']),
     'B': set(['C', 'G', 'T']),
     'V': set(['A', 'C', 'G']),
     'D': set(['A', 'G', 'T']),
     'N': set(['A', 'C', 'G', 'T']),
     '-': set(['A', 'C', 'G', 'T'])
 }


def count_wobbles(seq):
    return sum((1 for s in seq if s.upper() not in ['A', 'C', 'G', 'T']))


def hamming_dist(str1, str2):
    '''hamming distance'''
    assert len(str1) == len(str2), 'string lengths do not match'
    i = zip(str1, str2)
    n = 0
    while True:
        try:
            p = next(i)
        except StopIteration:
            break
        if p is None:
            break
        if (dna_code[p[0].upper()] & dna_code[p[1].upper()] == set([])):
            n += 1
    return n


def alignment_penalty(str1, str2):
    '''align with biopython'''
    str1, str2 = str1.upper(), str2.upper()
    score = pairwise2.align.globalxx(str1, str2, score_only = True)
    n = max(len(str1), len(str2)) - score
    return n


def trim_fwd_primer(seq, fwd, m):
    '''trim fwd primer'''
    for f in fwd:
        n = len(f)
        w = count_wobbles(f)
        if len(seq) >= n:
            #if hamming_dist(seq[:n], f) <= m:
            if alignment_penalty(seq[:n], f) <= m + w:
                seq = seq[n:]
                break
    return seq


def trim_rev_primer(seq, rev, m):
    '''trim rev primer'''
    for r in rev:
        n = len(r)
        w = count_wobbles(r)
        if len(seq) >= n:
            #if hamming_dist(seq[-n:], r) <= m:
            if alignment_penalty(seq[-n:], r) <= m + w:
                seq = seq[:-n]
                break
    return seq


def main():
    ### parameters
    random_nuleotides = 4
    allowed_mismatches = 3
    minimal_trimmed_nucleotides = 30
    maximal_trimmed_nucleotides = 60
    
    ### input files
    sequences, fwd, rev = SeqIO.parse(open(sys.argv[1]),'fasta'), SeqIO.parse(open(sys.argv[2]),'fasta'), SeqIO.parse(open(sys.argv[3]),'fasta')
    
    ### convert primers to str and reverse compelement reverse primer
    fwd, rev = [str(s.seq) for s in fwd], [str(s.seq.reverse_complement()) for s in rev]
    
    ### remove primers
    for fasta in sequences:
        name, sequence = fasta.id, str(fasta.seq)
        length = len(sequence)
        
        ### trim random nucleotides at 5' end of sequence
        sequence = sequence[random_nuleotides:]
        
        ### remove fwd and rev primers
        sequence = trim_fwd_primer(sequence, fwd, allowed_mismatches)
        sequence = trim_rev_primer(sequence, rev, allowed_mismatches)
        
        ### filter output
        if minimal_trimmed_nucleotides <= length - len(sequence) <= maximal_trimmed_nucleotides:
            print(">%s|trimmed_%s" % (name, length - len(sequence)))
            print(sequence)


if __name__ == '__main__':
    main()