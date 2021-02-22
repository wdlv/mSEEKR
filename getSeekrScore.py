import argparse
from seekr.fasta_reader import Reader
import corefunctions
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
        # remove N and make sure sequences only have AGCT
        seq=seq.replace("N","")
        seq_set = sorted(list(set(seq + "ATCG")))
        assert seq_set == ['A', 'C', 'G', 'T'], 'input sequences have bases not belong to ACGTN'

        # scale kmer count number to counts/kb of current sequence
        seq_length = len(seq)
        scaled_increment = 1000 / (seq_length - k + 1)
        for i in range(0,seq_length-k+1):
            kmerDict[seq[i:i+k]] += scaled_increment
        
        # add scaled 1 to kmer whose count is 0
        onerow = list(kmerDict.values())
        onerow = [scaled_increment if x==0 else x for x in onerow]
        # assign list value to numpy matrix's each row
        kmerDataMatrix[index] = np.asarray(onerow, dtype=np.float32)
    
    # log2 transform the count matrix
    kmerDataMatrix = np.log2(kmerDataMatrix)

    return kmerDataMatrix

# get a pearson correlation score of zscores calcuated from query and interesting fastas
def getPearsonCorrelation(oneSeq, queryMean, queryStd, interestingMean, interestingStd):
    queryZscore = (oneSeq - queryMean) / queryStd
    interestingZscore = (oneSeq - interestingMean) / interestingStd
    oneSeqPearCorr = stats.pearsonr(queryZscore, interestingZscore)[1]

    return oneSeqPearCorr



#Load arguments, see help= for explanation
parser = argparse.ArgumentParser()
parser.add_argument('--queryFasta', type=str,help='Path to query fasta file', required=True)
#parser.add_argument('--nullFasta', type=str,help='Path to null fasta file', required=True)
parser.add_argument('--interestingFasta', type=str,help='Path to interesting fasta file', required=True)
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

    #SEEKR fasta reader module
    # get mean and std; axis = 0 -> column | axis = 1 -> row

    # query fasta
    F = Reader(args.queryFasta)
    seqs = F.get_seqs()
    kmerDataMatrix = getSeqsKmerProcessedCounts(seqs, kVals, alphabet)
    queryMean = np.mean(kmerDataMatrix, axis = 0)
    queryStd = np.std(kmerDataMatrix, axis = 0)

    # interesting fasta
    F = Reader(args.interestingFasta)
    seqs = F.get_seqs()
    kmerDataMatrix = getSeqsKmerProcessedCounts(seqs, kVals, alphabet)
    interestingMean = np.mean(kmerDataMatrix, axis = 0)
    interestingStd = np.std(kmerDataMatrix, axis = 0)

    #print(queryMean, queryStd, interestingMean, interestingStd)

    # read in mSEEKR output dataframe
    mseekrdf = pd.read_csv(args.mSEEKRdataframeDir, sep="\t")

    hitsSeqs = list(mseekrdf['Sequence'])

    pearsonlist = []
    # iterate each hits and get each's pearson correlation score
    for singleSeq in hitsSeqs:
        oneSeq = getSeqsKmerProcessedCounts(singleSeq, kVals, alphabet)
        print(singleSeq)
        print(oneSeq)
        oneSeq = oneSeq[0]
        break
'''
        pearsonlist.append(getPearsonCorrelation(oneSeq, queryMean, queryStd, interestingMean, interestingStd))
    
    # add pearsonC list to mSEEKR dataframe
    mseekrdf['seekrPearsonCorr'] = pearsonlist

    mseekrdf.to_csv(f'./{args.name}_{kVals}_seekr_pearson_score.txt',sep='\t', index=False)

'''


