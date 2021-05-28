import argparse
import corefunctions
from scipy import stats
from itertools import product
import numpy as np
import pandas as pd
import gc


'''
Process background fasta file and get each sequence's kmer matrix then calcualte the pearson score between one sequence and anyone sequence in the fasta file
'''

#Load arguments, see help= for explanation
parser = argparse.ArgumentParser()
# parser.add_argument('--queryFasta', type=str,help='Path to query fasta file', required=True)
parser.add_argument('--backgroundFasta', type=str,help='Path to lncRNA background sequences fasta file; used to calculate mean and standard deviation of k-mer seekr score matrix', required=True)
# parser.add_argument('--mSEEKRdataframeDir',type=str,help='Directory to read in the output dataframe generated from command python mSEEKR.py',default='./', required=True)
parser.add_argument('--name',type=str,help='Name for output background fasta kmer seekr score matrix mean and std',default='bgMatrixMeanStd')
parser.add_argument('--outdir',type=str,help='Directory to save output background fasta kmer seekr score matrix mean and std',default='./')
parser.add_argument('-k',type=str,help='The same k value used in the python mSEEKR.py step', required=True)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG)',default='ATCG')

args = parser.parse_args()



if __name__ == '__main__':

    # Read in specified values of k, and the alphabet
    kVals = int(args.k)
    alphabet = args.a.upper()

    ### get mean and std; axis = 0 -> column | axis = 1 -> row

    # background fasta
    seqs = corefunctions.getCookedFasta(args.backgroundFasta)[1::2]
    backgroundKmerDataMatrix = corefunctions.getSeqsKmerProcessedCounts(seqs, kVals, alphabet)

    del seqs
    gc.collect()

    backgroundMean = np.mean(backgroundKmerDataMatrix, axis = 0)
    backgroundStd = np.std(backgroundKmerDataMatrix, axis = 0)

    backgroundKmerDataMatrix = (backgroundKmerDataMatrix - backgroundMean) / backgroundStd

    backgroundseekrScoreMatrix = corefunctions.getSeekrScorePearson(backgroundKmerDataMatrix, backgroundKmerDataMatrix)

    del backgroundKmerDataMatrix
    gc.collect()

    # get only the upper triangle of the matrix
    target_values = backgroundseekrScoreMatrix[np.triu_indices_from(backgroundseekrScoreMatrix, k=1)]

    del backgroundseekrScoreMatrix
    gc.collect()

    target_mean = np.mean(target_values)
    target_std = np.std(target_values)

    output = (target_mean, target_std)

    np.savetxt(args.outdir + "/" + args.name + ".txt", output)



