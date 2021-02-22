import argparse
from seekr.fasta_reader import Reader
#import corefunctions
from scipy import stats
from itertools import product
import numpy as np
import pandas as pd

# preprocess the kmer counts
def getSeqsKmerProcessedCounts(seqs, k, alphabet):
    # get all the possible kmer combination
    demoKmers = [''.join(p) for p in product(alphabet, repeat=k)]

    # create a numpy array with initial value 0; 4 is equal to len(ATCG)
    kmerDataMatrix = np.zeros((len(seqs), 4**k), dtype=np.float32)

    for index, seq in enumerate(seqs):
        # create a dict with all kmers and initial value is 0
        kmerDict = dict.fromkeys(demoKmers, 0)

        # scale kmer count number to counts/kb of current sequence
        seq_length = len(seq)
        scaled_increment = 1000 / (seq_length - k + 1)
        for i in range(0,seq_length-k+1):
            currentKmer = seq[i:i+k]
            if currentKmer in kmerDict:
                kmerDict[currentKmer] += scaled_increment
        
        # add scaled 1 to kmer whose count is 0
        onerow = list(kmerDict.values())
        onerow = [scaled_increment if x==0 else x for x in onerow]
        # assign list value to numpy matrix's each row
        kmerDataMatrix[index] = np.asarray(onerow, dtype=np.float32)
    
    # log2 transform the count matrix
    kmerDataMatrix = np.log2(kmerDataMatrix)

    return kmerDataMatrix

# get a pearson correlation score of zscores calcuated from query and background fastas
def getPearsonCorrelation(oneSeq, querySeq, backgroundMean, backgroundStd):
    queryZscore = (querySeq - backgroundMean) / backgroundStd
    backgroundZscore = (oneSeq - backgroundMean) / backgroundStd
    oneSeqPearCorr = stats.pearsonr(queryZscore, backgroundZscore)[1]

    return oneSeqPearCorr



#Load arguments, see help= for explanation
parser = argparse.ArgumentParser()
parser.add_argument('--queryFasta', type=str,help='Path to query fasta file', required=True)
#parser.add_argument('--nullFasta', type=str,help='Path to null fasta file', required=True)
parser.add_argument('--backgroundFasta', type=str,help='Path to lncRNA background sequences fasta file; used to calculate mean and standard deviation for each k-mer', required=True)
parser.add_argument('--mSEEKRdataframeDir',type=str,help='directory to read in mSEEKR dataframe',default='./', required=True)
parser.add_argument('--name',type=str,help='name for seekPearsonCorr file',default='out', required=True)
parser.add_argument('--dir',type=str,help='directory to save output',default='./')
parser.add_argument('-k',type=str,help='Comma delimited string of possible k-mer values', required=True)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')

args = parser.parse_args()


if __name__ == '__main__':

    # Read in specified values of k, and the alphabet
    kVals = int(args.k)
    alphabet = args.a.upper()

    # SEEKR fasta reader module
    # get mean and std; axis = 0 -> column | axis = 1 -> row

    # query fasta. Assume query fasta is one sequence. If not, merge multiple sequences to one sequence
    F = Reader(args.queryFasta)
    seqs = F.get_seqs()

    if len(seqs) > 1:
        seqs = ''.join(seqs)
        querySeq = getSeqsKmerProcessedCounts([seqs], kVals, alphabet)[0]
    querySeq = getSeqsKmerProcessedCounts(seqs, kVals, alphabet)[0]

    # background fasta
    F = Reader(args.backgroundFasta)
    seqs = F.get_seqs()
    kmerDataMatrix = getSeqsKmerProcessedCounts(seqs, kVals, alphabet)
    backgroundMean = np.mean(kmerDataMatrix, axis = 0)
    backgroundStd = np.std(kmerDataMatrix, axis = 0)

    # read in mSEEKR output dataframe
    mseekrdf = pd.read_csv(args.mSEEKRdataframeDir, sep="\t")

    hitsSeqs = list(mseekrdf['Sequence'])

    pearsonlist = []
    # iterate each hits and get each's pearson correlation score
    for singleSeq in hitsSeqs:
        oneSeq = getSeqsKmerProcessedCounts([singleSeq], kVals, alphabet)[0]
        pearsonlist.append(getPearsonCorrelation(oneSeq, querySeq, backgroundMean, backgroundStd))
    
    # add pearsonC list to mSEEKR dataframe
    mseekrdf['seekrPearsonCorr'] = pearsonlist

    mseekrdf.to_csv(f'./{args.name}_{kVals}_seekr_pearson_score.txt',sep='\t', index=False)

