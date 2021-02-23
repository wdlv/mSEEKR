import corefunctions
#import coreStats
import argparse
import itertools
import numpy as np
from seekr.fasta_reader import Reader
from scipy.stats import norm
from collections import defaultdict
#from multiprocessing import pool
from scipy.stats import gaussian_kde
from itertools import product
import sys
import os
import glob
import pickle

'''
This program calcuate the emission matrix based on kmer counts. Prepare transition matrix, states and 
saves them to a binary file that countains a dictionary of dictionaries
------------------------------------------------------------------------
query: str
    query sequence. Represented by '+' in certain code
null: str
    null sequence. Represented by '-' in certain code
qT: float
    transition rate from query to query so transition rate from query to null would be 1 - qT
nT: float
    transition rate from null to null so transition rate from null to query would be 1 - nT
alphabet: list
    a list of base pairs. for example ['A', 'T', 'C', 'G']
qCount: dict
    query sequence reads' kmer count. 
    data structure example: 
    {2: {'AA': 30, 'AT': 24, 'AC': 18,.....}, 3: {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....}, 4: {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}}
nCount: dict
    null sequence reads' kmer count.
    data structure example:
    {2: {'AA': 30, 'AT': 24, 'AC': 18,.....}, 3: {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....}, 4: {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}}
kVals: list
    a list of k values - for example [2, 3, 4]
qKCount: dict
    query sequence's specific k's kmer counts dictionary
    data structure example:
    {'AA': 30, 'AT': 24, 'AC': 18,.....} or {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....} or {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}
nKCount: dict
    null sequence's specific k's kmer counts dictionary
    data structure example:
    {'AA': 30, 'AT': 24, 'AC': 18,.....} or {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....} or {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}
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
kmers: list
    a list of specific k length kmers
    Format example:
    k=3 ['AAA', 'AAT', 'AAC', 'AAG', 'ATA'....]
    k=4 ['AAAA', 'AAAT', 'AAAC', 'AAAG', 'AATA',....]
'''

# Initialize program arguments, see help= for explanation of each
parser = argparse.ArgumentParser()
parser.add_argument("--query",type=str,help='Path to kmer count file for sequences of interest (e.g. functional regions of a ncRNA)', required=True)
parser.add_argument('--null', type=str,help='Path to kmer count file that compose null model (e.g. transcriptome, genome, etc.)', required=True)
parser.add_argument('--qT',type=float,help='Probability of query to query transition', required=True)  #default=.999)
parser.add_argument('--nT',type=float,help='Probability of null to null transition', required=True) #default=.9999)
parser.add_argument('--qPrefix',type=str,help='String, Output file prefix;default=None',default='query', required=True)
parser.add_argument('--nPrefix',type=str,help='String, Output file prefix;default=None',default='null', required=True)
parser.add_argument('--dir',type=str,help='Output directory',default='./', required=True)
parser.add_argument('-k',type=str,help='Comma delimited string of possible k-mer values,must be found in the k-mer count file', required=True) #default='2,3,4')
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')

# get input argument
args = parser.parse_args()
kVals = [int(i) for i in args.k.split(',')]

# Check if specified directory exists
# If yes, prompt if replace, else end program
# If no, crash
# Else, loop
# create new directory if not existing
if not args.dir.endswith('/'):
    args.dir+='/'
newDir = f'{args.dir}{args.qPrefix}_{args.nPrefix}/'

if not os.path.exists(newDir):
    os.mkdir(newDir)
else:
    flag = True
    while flag:
        usrIN = input(f'Directory {newDir} exists, continue? y/n: ').strip().lower()
        if usrIN == 'y':
            flag = False
        elif usrIN == 'n':
            print('Initiating self-destruct sequence')
            sys.exit()
        else:
            print('Please enter y or n')


alphabet = list(args.a)  # like ['A', 'T', 'C', 'G']

# Load k-mer counts
qCount = pickle.load(open(args.query,'rb'))
nCount = pickle.load(open(args.null,'rb'))


# Loop through specified values of k
# Check if they exist in the counts file,
# and call corefunctions.HMM to generate the HMM matrices
for k in kVals:
    if (k in qCount.keys()) and (k in nCount.keys()):
        qKCount = qCount[k]
        nKCount = nCount[k]
        kDir = newDir+f'{k}/'
        if not os.path.exists(kDir):
            os.mkdir(kDir)
        A,E,states,pi = corefunctions.HMM(qKCount,nKCount,k,args.a,args.qT,args.nT)
        # kmers = [''.join(p) for p in itertools.product(alphabet,repeat=k)]
        # queryMkv = corefunctions.transitionMatrix(qKCount,k,alphabet)
        # nullMkv = corefunctions.transitionMatrix(nKCount,k,alphabet)
        # lgTbl = corefunctions.logLTbl(queryMkv,nullMkv)
    else:
        print(f'Missing {k}-mer counts in count file... skipping')

    # np.savetxt(f'{kDir}logtbl.mkv',lgTbl)
    pickle.dump({'A':A,'E':E,'pi':pi,'states':states},open(f'{kDir}hmm.dict','wb'))
