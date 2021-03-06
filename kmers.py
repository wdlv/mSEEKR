import kmers
import argparse
import pickle
import os
# from multiprocessing import pool
from itertools import starmap
from itertools import product
import corefunctions

'''
This program counts k-mers for multiple specified values of k and saves them
to a binary file that countains a dictionary of dictionaries
------------------------------------------------------------
args: object 
    contains all the input options
kVals: list
    a lot of different integral k values
a: str
    a short string to generate kmers. default is ATCG
F: Reader object
    an object of Reader class from package seekr
fS: list
    a list of sequence reads without header
    Format example:
    ['AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCC', 'CTCCGCGTGGTCTATGATGGTGCATTTTGGTCCAGTCAGGCCCGGTGTGG', 'TCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGGCGTGGTTACG',....]
fString: str
    a string of sequence reads with $ being the delimiter character
lenFString: int
    number of total base pairs
dataDict: dict
    data structure example: 
    {2: {'AA': 30, 'AT': 24, 'AC': 18,.....}, 3: {'AAA': 23, 'AAT': 2, 'AAC': 7, 'AAG': 1....}, 4: {'AAAA': 19, 'AAAT': 2, 'AAAC': 4,....}}
kDir: str
    path to store generated output file, input by the user

'''

#Load arguments, see help= for explanation
parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=str,help='Path to fasta file', required=True)
parser.add_argument('--name',type=str,help='Desired output name for count file',default='out', required=True)
parser.add_argument('--dir',type=str,help='Directory to save output count file',default='./')
parser.add_argument('-k',type=str,help='Comma delimited string of possible k-mer values. For example, 3,4,5 or just 4', required=True)
parser.add_argument('-a',type=str,help='String, Alphabet to generate k-mers (e.g. ATCG); default=ATCG',default='ATCG')
# parser.add_argument('-n',type=int,help='Number of CPU cores. Each job corresponds to a value of k, and the program scales well with multiprocessing',default=1)

args = parser.parse_args()



if __name__ == '__main__':

    # Read in specified values of k, and the alphabet
    kVals = [int(i) for i in args.k.split(',')]
    a = args.a.upper()

    #Read in fasta file
    seqs = corefunctions.getCookedFasta(args.fasta)[1::2]


    #Join sequences together using $ delimiter character
    fString = '$'.join(seqs)
    #lenFString = sum([len(i) for i in fS])

    # Need to figure out how to deal with very long fasta files (~ 2-3X the size of the transcriptome in mice)
    # if lenFString >= 2147483647:
    #     fString='$'.join(fS[::10]).upper()

    #Split jobs onto processors and call kmers.pyx cython file
    # with pool.Pool(args.n) as multiN:
    #     jobs = multiN.starmap(kmers.main,product(*[[fString],kVals,[a]]))
    #     dataDict = dict(jobs)
    #     print(dataDict.keys())

    # call kmers.pyx cython file to conduct kmer calculation and get kmer count dictionary as return
    dataDict = dict(starmap(kmers.main,product(*[[fString],kVals,[a]])))


    #Save data
    kDir = args.dir
    if not kDir.endswith('/'):
        kDir+='/'
    pickle.dump(dataDict,open(f'{kDir}{args.name}.dict','wb'))

