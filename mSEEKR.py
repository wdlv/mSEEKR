import corefunctions
import argparse
from itertools import product
import numpy as np
from collections import defaultdict
from itertools import starmap
#from multiprocessing import pool
import sys
import os
import pickle
from math import log
import pandas as pd
from operator import itemgetter

'''
This program use precalculated emission matrix, prepared transition matrix, pi and states to find out HMM state path through sequences of interest
-------------------------------------------------------------------------------------------------------------------------------------------------
alphabet: list
    a list of base pairs. for example ['A', 'T', 'C', 'G']
model: str
    Path to .mkv file output from train.py or bw.py
A: dict
    Hidden state transition matrix
    initial value is: 
    {'+':{'+':np.log2(args.qT),'-':np.log2(1-args.qT)},'-':{'+':np.log2(1-args.nT),'-':np.log2(args.nT)}}
E: dict
    Hidden state emission matrix
    Format example:
    {'+': {'AAAA': -5.082989364671032, 'AAAT': -8.330916878114618, 'AAAC': -7.330916878114617,.....'GGGC': -6.523561956057013, 'GGGG': -6.523561956057013}, '-': {'AAAA': -6.735347642028761, 'AAAT': -7.242499465916189, 'AAAC': ...}}
states: dict
    states = ('+','-')
pi: dict
    Starting probability of each hidden state (+ or -)
    initial value is: 
    {'+':np.log2(.5),'-':np.log2(.5)} 
hmm: dict
    Format example:
    {
    'A':{'+': {'+': -0.00014427671804501932, '-': -13.287712379549609}, '-': {'+': -13.287712379549609, '-': -0.00014427671804501932}},
    'E':{'+': {'AAAA': -5.082989364671032, 'AAAT': -8.330916878114618, 'AAAC': -7.330916878114617, ....}, '-': {'AAAA': -6.735347642028761, 'AAAT': -7.242499465916189, 'AAAC': ...}}
    'pi':{'+': -1.0, '-': -1.0},
    'states':('+', '-')
    }
O: list
    observed sequence of k-mers. (doesn't have 'N' in kmers)
kmers: list
    a list of specific k length kmers
    Format example:
    k=3 ['AAA', 'AAT', 'AAC', 'AAG', 'ATA'....]
    k=4 ['AAAA', 'AAAT', 'AAAC', 'AAAG', 'AATA',....]
target: Reader object
    an object of Reader class from package seekr
targetSeqs: list
    a list of sequence reads without header
    Format example:
    ['AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCC', 'CTCCGCGTGGTCTATGATGGTGCATTTTGGTCCAGTCAGGCCCGGTGTGG', 'TCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGGCGTGGTTACG',....]
targetHeaders: list
    a list of headers
    Format example
    ['>ENSMUST00000193812.1|ENSMUSG00000102693.1|OTTMUSG00000049935.1|OTTMUST00000127109.1|RP23-271O17.1-001|RP23-271O17.1|1070|', '>ENSMUST00000195335.1|ENSMUSG00000103377.1|OTTMUSG00000049960.1|OTTMUST00000127145.1|RP23-317L18.4-001|RP23-317L18.4|2819|'......]
dataDict: dict
    key is a sequence header
    value is a dataframe
    Format example:
    {'>mm10kcnq1ot1':     Start    End    kmerLLR        seqName                                           Sequence
    0     835    999  84.173049  >mm10kcnq1ot1  ATTCGTGCCGCGCTTTCGCGGCTGGGCTCCATCTTCGTTTTGCCGC...
    1   79184  79257  72.086333  >mm10kcnq1ot1  ATTATTTTGTGTCTTTTTTTGTTTGTTTGTTTTTTGTTTTTTGTTT...
    2   69556  69605  63.331314  >mm10kcnq1ot1  TCCCAACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA.....}
O: list
    a list of k length kmers extracted from sequence and exclude all kmers with 'N' inside.
    Format example:
    when k=4, O is something like ['GGAC', 'GACA', 'ACAG', 'CAGC', 'AGCA',...]
oIdx: list
    a list of index. These index is the real location of O's kmers in the original sequence
    Format example:
    [0, 1, 2, 3, 6, 7, 8, 9...]
nBP: list
    a list of tuples.
    Format example:
    [(36, 'N'), (111, 'N'), (300, 'N').....]
bTrack: list
    a list of query/null states
    Format example:
    [-', '-', '-', '-', '-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',....]

'''


''' hmmCalc
Run several functions including viterbi algorithm, log-likelihood, and generate output dataframes
Input: fasta file information
Output: dataframe object

'''
def hmmCalc(data):
    tHead,tSeq = data
    O,oIdx,nBP = corefunctions.kmersWithAmbigIndex(tSeq,k)
    A,E,states,pi= hmm['A'],hmm['E'],hmm['states'],hmm['pi']
    bTrack = corefunctions.viterbi(O,A,E,states,pi)
    #Zip the indices of unambig k-mers with their viterbi derived HMM state labels
    coordBTrack = list(zip(oIdx,bTrack)) # [(1,'-'),(2,'+',...(n,'+'))]
    mergedTrack = coordBTrack + nBP # add back in ambig locations
    mergedTrack.sort(key=itemgetter(0)) # sort master list by index
    hmmTrack = [i[1] for i in mergedTrack] # fetch just state label from mergedTrack ['-','+',...,'+']
    groupedHits = corefunctions.groupHMM(hmmTrack) # ['-----','++++++++++','-','++++','------------']

    # Return sequences of HMM hits, and their start and end locations in the original sequence
    seqHits,starts,ends = corefunctions.formatHits(groupedHits,k,tSeq)
    if (seqHits):
        df = corefunctions.hitOutput(seqHits,starts,ends,k,E,tHead,tSeq)
        return tHead,df
    # Alternative output (transcript by transcript)

    else:
        return tHead,None


# input commands
parser = argparse.ArgumentParser()
parser.add_argument("--model",type=str,help='Path to .dict file output from train.py or bw.py', required=True)
parser.add_argument("-k",type=int,help='Value of k to use. Must be the same as the k value used in training', required=True)
parser.add_argument('--db',type=str,help='Path to fasta file with sequences to calculate similarity score', required=True)
parser.add_argument('--prefix',type=str,help='String, Output file name', required=True)
#parser.add_argument('--bkg',type=str,help='Path to fasta file from which to calculate background nucleotide frequencies, if not passed default is uniform',default=None)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')
#parser.add_argument('-n',type=int,help='Integer 1 <= n <= max(cores), Number of processor cores to use; default = 1. This scales with the number of sequence comparisons in --db',default=1)
parser.add_argument('--fasta',action='store_true',help='FLAG: print sequence of hit, ignored if --wt is passed')

args = parser.parse_args()


if __name__ == '__main__':

    #Loop over values of k
    args.a = args.a.upper()
    model = args.model

    hmm = pickle.load(open(model,'rb'))
    A,E,pi,states = hmm['A'],hmm['E'],hmm['pi'],hmm['states']

    # Explicitly determine k from the size of the log matrix and the size of the alphabet used to generate it
    k = int(log(len(hmm['E']['+'].keys()),len(args.a)))
    assert k == args.k, 'Value of k provided does not match supplied hmm file'

    cookedFasta = corefunctions.getCookedFasta(args.db)
    targetSeqs,targetHeaders = cookedFasta[1::2],cookedFasta[::2]

    #Pool processes onto number of CPU cores specified by the user
    # with pool.Pool(args.n) as multiN:
    #     jobs = multiN.starmap(hmmCalc,product(*[list(zip(targetHeaders,targetSeqs))]))
    #     dataDict = dict(jobs)

    dataDict = dict(starmap(hmmCalc,product(*[list(zip(targetHeaders,targetSeqs))])))


    #Check if no hits were found
    # if not all(v == None for v in dataDict.values()):
    dataFrames = pd.concat([df for df in dataDict.values() if not None])
    dataFrames['Start']+=1 #1-start coordinates
    dataFrames['End']
    dataFrames['Length'] = dataFrames['End'] - dataFrames['Start'] +1
    dataFrames = dataFrames[['Start','End','Length','kmerLLR','seqName','Sequence']]
    if not args.fasta:
        dataFrames = dataFrames[['Start','End','Length','kmerLLR','seqName']]
    dataFrames.sort_values(by='kmerLLR',ascending=False,inplace=True)
    dataFrames.reset_index(inplace=True,drop=True)
    dataFrames.to_csv(f'./{args.prefix}_{k}_viterbi.txt',sep='\t', index=False)
