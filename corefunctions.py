import numpy as np
from itertools import product
from itertools import groupby
from scipy.special import logsumexp
import pandas as pd

'''
transform sequences to kmer counts
'''
# preprocess the kmer counts
def getSeqsKmerProcessedCounts(seqs, k, alphabet):
    # get all the possible kmer combination
    demoKmers = [''.join(p) for p in product(alphabet, repeat=k)]

    # create a numpy array with initial value 0; 4 is equal to len(ATCG)
    kmerDataMatrix = np.zeros((len(seqs), 4**k), dtype=np.float32)

    for index, seq in enumerate(seqs):
        # create a dict with all kmers and initial value is 1
        kmerDict = dict.fromkeys(demoKmers, 1)

        # scale kmer count number to counts/kb of current sequence
        seq_length = len(seq)
        scaled_increment = 1000 / (seq_length - k + 1)
        #scaled_increment = 1 if seq_length <= 1000 else scaled_increment

        for i in range(0,seq_length-k+1):
            currentKmer = seq[i:i+k]
            if currentKmer in kmerDict:
                kmerDict[currentKmer] += scaled_increment
        
        onerow = list(kmerDict.values())


        # add scaled 1 to kmer whose count is 0
        # onerow = [scaled_increment if x==0 else x for x in onerow]

        # assign list value to numpy matrix's each row
        kmerDataMatrix[index] = np.asarray(onerow, dtype=np.float32)
    
    # log2 transform the count matrix
    kmerDataMatrix = np.log2(kmerDataMatrix)

    return kmerDataMatrix

# get a pearson correlation score of zscores calcuated from query and background fastas
# def getPearsonCorrelation(oneSeq, querySeq, backgroundMean, backgroundStd):

#     queryZscore = (querySeq - backgroundMean) / backgroundStd
#     backgroundZscore = (oneSeq - backgroundMean) / backgroundStd
#     oneSeqPearCorr = stats.pearsonr(queryZscore, backgroundZscore)[0]

#     return oneSeqPearCorr

'''
The following function getSeekrScorePearson is using the same strategy as the method used in seekr package
'''

# def rowNormalization(inputSeqs):

#     inputSeqs = (inputSeqs.T - np.mean(inputSeqs, axis=1)).T
#     inputSeqs = (inputSeqs.T / np.std(inputSeqs, axis=1)).T

#     return inputSeqs

def getSeekrScorePearson(inputSeqs1, inputSeqs2):
    # conduct row normalization to get normalized matrix
    inputSeqs1 = (inputSeqs1.T - np.mean(inputSeqs1, axis=1)).T
    inputSeqs1 = (inputSeqs1.T / np.std(inputSeqs1, axis=1)).T
    inputSeqs2 = (inputSeqs2.T - np.mean(inputSeqs2, axis=1)).T
    inputSeqs2 = (inputSeqs2.T / np.std(inputSeqs2, axis=1)).T
    
    # seekr's pearson
    return np.inner(inputSeqs1, inputSeqs2) / inputSeqs1.shape[1]


'''
Read in fasta file and reorganize it to a list with following format:
[header1, sequence1, header2, sequence2, .....]
'''

def getCookedFasta(fastaFile):
    with open(fastaFile) as ff:
        rawFasta = [line.strip() for line in ff]
    rawFasta = [x for x in rawFasta if len(x)>0]
    cookedFasta = []
    seq = ""
    for i, line in enumerate(rawFasta):
        if line[0] == ">":
            if seq:
                cookedFasta.append(seq.upper())
                seq = ""
            else:
                assert i == 0, "There may be a header without a sequence at line {}.".format(i)
            cookedFasta.append(line)
        else:
            seq += line
    cookedFasta.append(seq.upper())
    return cookedFasta



''' Key for itertools groupby
    Alters flag when sequence changes from one condition to another
    Input: Sequence of characters with some alphabet and a trigger condition
    Output: flag: 0 or 1
'''

class Key(object):
    def __init__(self):
        self.is_nt,self.flag,self.prev = ['-','N'],[0,1],None
    def __call__(self,e):
        # Initial True/False  if first char in string is + or -
        ebool = any(x in self.is_nt for x in e)
        # If key exists (self.prev is defined), do true/false check
        # else, set value to false
        if self.prev:
            prevbool = any(x in self.is_nt for x in self.prev)
        else:
            prevbool = None
        # if string goes from - to +, or + to -, swap flag
        if prevbool != ebool:
            self.flag = self.flag[::-1]
        # set previous encountered char, for the next interation of this, to the current value
        self.prev = e
        return self.flag[0]




''' groupHMM
Return a list of strings separating HMM state labels
Input: String
Output: List of lists
Output example: ['---','++','------','+','-------',...]
'''
def groupHMM(seq):
    return [''.join(list(g)) for k,g in groupby(seq,key=Key())] # Key() refer to class Key(object) above





'''
Combine all the data to create a dataframe
E: emission matrix
seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
'''

def hitOutput(seqHits,starts,ends,k,E,tHead,tSeq):
    info = list(zip(seqHits,starts,ends)) # example [('GGCCCGGTGTGGTCGGCCTCATTTTGGAT.......', 88, 177),......]
    dataDict = dict(zip(list(range(len(seqHits))),info))
    df = pd.DataFrame.from_dict(dataDict,orient='index')
    #calculate log-likelihood ratio of k-mers in the + model vs - model
    df['kmerLLR'] = LLR(seqHits,k,E)
    df['seqName'] = tHead
    df.columns = ['Sequence','Start','End','kmerLLR','seqName']
    df.sort_values(by='kmerLLR',inplace=True,ascending=False)
    df.reset_index(inplace=True)
    fa = df['Sequence']
    df = df[['Start','End','kmerLLR','seqName','Sequence']]

    return df



''' LLR
Return log-likelihood ratio between two models in HMM for + k-mers
Input: sequnce of hits, value of k, k-mer frequencies in HMM emmission matrix
Output: Array of LLRs for each hit

hits = seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
'''
def LLR(hits,k,E):
    arr = np.zeros(len(hits))
    for i,hit in enumerate(hits):
        LLRPos,LLRNeg=0,0
        for j in range(len(hit)-k+1):
            kmer=hit[j:j+k]
            LLRPos += E['+'][kmer]
            LLRNeg += E['-'][kmer]
        llr = LLRPos-LLRNeg
        arr[i] = llr
    return arr


# def getRawKmerCount(seq, k):
#     kmerCountDict = defaultdict(int)
#     seq_length = len(seq)
#     for i in range(0,seq_length-k+1):
#         currentKmer = seq[i:i+k]
#         if 'N' not in currentKmer:
#             kmerCountDict[currentKmer] += 1
#     print(kmerCountDict)

# def transcriptOutput(seqHits,starts,ends,k,E,tHead,tSeq):

#     sumHits = LLR(seqHits,k,E)
#     lens = ends-starts # lengths of individual hits
#     df = pd.DataFrame([np.sum(sumHits)]) # sum of hits
#     df['totalLenHits'] = (np.sum(lens)) # sum of all hit lengths
#     df['fracTranscriptHit'] = df['totalLenHits']/len(tSeq) # fraction of transcript that is hit
#     df['longestHit'] = np.max(lens) # longest HMM hit
#     df['seqName'] = tHead
#     df.columns = ['sumLLR','totalLenHits','fracTranscriptHit','longestHit','seqName']

#     return df




''' kmersWithAmbigIndex
Return list of kmers, indices of k-mers without ambiguity, and indices of those
with ambiguity
Input: string
Output: List of string, list of indices, list of indices
'''
def kmersWithAmbigIndex(tSeq,k):
    O = [tSeq[i:i+k].upper() for i in range(0,len(tSeq)-k+1)]
    O = [o for o in O if 'N' not in o]  # example: when k=4, O is something like ['GGAC', 'GACA', 'ACAG', 'CAGC', 'AGCA',...]
    # Match k-mers without ambig char to index in original string
    oIdx = [i for i in range(0,len(tSeq)-k+1) if 'N' not in tSeq[i:i+k]] # example: when k=4, oIdx is something like [0, 1, 2, 3, 4, 5, 6, 7, 8, 9...]
    # Match k-mers with ambig char to index in original string
    nBP = [i for i in range(0,len(tSeq)-k+1) if 'N' in tSeq[i:i+k]]
    # zip the indices with marker character N
    nBP = list(zip(nBP,['N']*len(nBP))) # example. nBP is something like [(36, 'N'), (111, 'N'), (300, 'N').....]
    return O, oIdx, nBP





# def getFwd(seqHits,A,pi,states,E,k,alphabet):
#     fwdPs = []
#     for hit in seqHits:
#         O = [hit[i:i+k].upper() for i in range(0,len(hit)-k+1)]
#         O = [o for o in O if 'N' not in o]
#         '''
#         forward algorithm to calculate log P(O|HMM)
#         '''
#         fP = fwd(O,A,pi,states,E,k,alphabet)
#         fwdPs.append(fP)
#     return fwdPs





'''
Process raw results to find out the break point of the original sequence then break it to corresponding query reads or null reads fragments based on groupedHits. 
Also the start position and end position of each fragments in the original sequence are recorded.

From:
groupedHits: ['-----','++++++++++','-','++++','------------']
tSeq: 'AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCCCTCCGCGTGG......'

to:
seqHits: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
starts: [   88   498   835  5614 14919 15366 19473 22753 69402 69556.....]
ends: [  177   627   999  5636 14971 15388 19502 22776 69436 69605.......]
'''
def formatHits(groupedHits,k,tSeq):
    idx = 0
    indexGroupHits = []
    # Loop below formats the hmm output as such:
    # [([0,1,2]),'---'),([3,4],'++'),([5],'-'),...]
    # Grouping HMM states with their correct index in the list of k-mers
    for i,group in enumerate(groupedHits):   # groupedHits example ['-----','++++++++++','-','++++','------------']
        indexGroupHits.append([])
        for kmer in group:
            indexGroupHits[i].append(idx)
            idx+=1
    hits = list(zip(indexGroupHits,groupedHits)) # hits example [([0,1,2]),'---'),([3,4],'++'),([5],'-'),...]
    
    seqHits = []
    seqHitCoords = []
    for group in hits:
        if '+' in group[1]:
            start,end = group[0][0],group[0][-1]+k #convert k-mer coord to bp coord
            seqHitCoords.append(f'{start}:{end}')
            seqHits.append(tSeq[start:end])
    # seqHits - squence hits base pair version - example: ['GGCCCGGTGTGGTCGGCCTCATTTTGGATTACTTCGGTGGGCTTCTCCTCGG...', 'TCCGTGTGGATCGTTTCAGCACGGATC......'......]
    # tSeq example 'AGGAGGAACAGTTGCCTCAGCACGTCTGCGCAGCTTTCCTTGCGGCGCCCCTCCGCGTGG......'
    # seqHitCoords example ['88:177', '498:627', '835:999', '5614:5636', '14919:14971', '15366:15388', '19473:19502', '22753:22776', '69402:69436', '69556:69605'......]
    starts = np.array([int(c.split(':')[0]) for c in seqHitCoords]) # starts example [   88   498   835  5614 14919 15366 19473 22753 69402 69556.....]
    ends = np.array([int(c.split(':')[1]) for c in seqHitCoords]) # ends example [  177   627   999  5636 14971 15388 19502 22776 69436 69605.......]
    return seqHits,starts,ends





# def transitionMatrix(kmers,k,alphabet):
#     states = np.zeros((4**(int(k)-1), 4), dtype=np.float64)
#     stateKmers = [''.join(p) for p in product(alphabet,repeat=k-1)]
#     for i, currState in enumerate(stateKmers):
#         tot = 0
#         for j, nextState in enumerate(alphabet):
#             count = kmers[currState+nextState]
#             tot += count
#         if tot > 0:
#             for j, nextState in enumerate(alphabet):
#                 states[i, j] = kmers[currState+nextState] / float(tot)
#     return states

# # calculate nucleotide content of a sequence... unused now
# def nucContent(nullSeqs,alphabet):
#     seqs = ''.join(nullSeqs)
#     seqs = seqs.upper()
#     freq = [seqs.count(nt)/len(seqs) for nt in alphabet]
#     return dict(zip(alphabet,freq))

# # Calculate the log2 odds table between two matrices
# def logLTbl(q,null):
#     return np.log2(q) - np.log2(null)




'''
HMM: Generate dictionaries containing information necessary to perform the viterbi algorithm

Inputs: qKCounts - Same certain length's Kmers' query count dictionary
        nKCounts - Same certain length's Kmers' null count dictionary
        k - value of k
        alphabet - alphabet (ATCG)
        qT: + to + transition probability (query to query, query to null is 1-qT)
        nT: - to - transition probability (null to null, null to query is 1-nT)
Returns:    A - Dictionary, Hidden state transition matrix
            E - Dictionary, Emission matrix, kmer counts
            state - tuple, list of states (+,-)
            pi - Dictionary, initial probability of being in + or -
'''
def HMM(qKCounts,nKCounts,k,alphabet,qT,nT):
    kmers = [''.join(p) for p in product(alphabet,repeat=k)]
    hmmDict = {}
    countArr = np.array(list(qKCounts.values()))
    # Convert raw counts to frequencies, then log probability
    hmmDict['+'] = np.log2(countArr/np.sum(countArr))
    #hmmDict['+'] = hmmDict['+'] - np.mean(hmmDict['+'])
    countArr = np.array(list(nKCounts.values()))
    # Convert raw counts to frequencies, then log probability
    hmmDict['-'] = np.log2(countArr/np.sum(countArr))
    #hmmDict['-'] = hmmDict['-'] - np.mean(hmmDict['-'])

    states = ('+','-')
    pi = {'+':np.log2(.5),'-':np.log2(.5)}
    A = {'+':{'+':np.log2(qT),'-':np.log2(1-qT)},'-':{'+':np.log2(1-nT),'-':np.log2(nT)}}
    E = {'+': dict(zip(kmers,hmmDict['+'])),'-':dict(zip(kmers,hmmDict['-']))}
    return A,E,states,pi




'''
Viterbi: Calculate the most likely sequence of hidden states given observed sequence, transition matrix, and emission matrix

Inputs: O - list, observed sequence of k-mers
        A - dictionary, transition matrices of hidden states
        E - dictionary, emission matrices of hidden states
        states - list, hidden states (+,-)
        m: + to + transition probability
        n: - to - transition probability
        pi - Dictionary, initial probability of being in + or -
Returns:    backTrack - list, sequence of hidden states

'''
def viterbi(O,A,E,states,pi):

    # Initialize list of dictionaries for the current step
    # and ukprev, which tracks the state transitions that maximize the 'viterbi function'
    uk=[{}]
    ukprev = [{}]
    N = len(O)
    # calculate initial probabilities in each state given the first kmer
    for state in states:
        uk[0][state]=pi[state]+E[state][O[0]]  # uk[0]['+'] = np.log2(.5) +  E['+']['AAAA']
        ukprev[0][state] = None # previous state does not exist, set to None
    # Loop through observed sequence
    # For each state, calculate the cumulative probability recursively
    # Store the state transition that maximizes this probability in ukprev for each current state
    for n in range(1,N):
        uk.append({})
        ukprev.append({})
        for state in states:
            prevSelState = states[0] # this is just an arbitrary choice to start checking at the start of the list
            currMaxProb = A[state][prevSelState] + uk[n-1][prevSelState] # probability function
            for pState in states[1:]: # now check the other states...
                currProb = A[state][pState] + uk[n-1][pState]
                if currProb > currMaxProb: # if larger then the first one we checked, we have a new winner, store and continue loop and repeat
                    currMaxProb = currProb
                    prevSelState = pState
            # The emission probability is constant so add at the end rather than in the loop
            max_prob = currMaxProb + E[state][O[n]]
            # save the cumalitive probability for each state
            uk[n][state] = max_prob
            # save the previous state that maximized this probability above
            ukprev[n][state] = prevSelState

    z = max(uk[-1],key=uk[-n].get) # retrieve the state associated with largest log probability
    prev = ukprev[-1][z] # get the state BEFORE "z" above that maximized z
    backtrack = [z,prev] # start our backtrack with knowns
    # Loop through BACKWARDS, getting the previous state that yielded the 'current' max state
    for n in range(N-2,-1,-1):
        backtrack.append(ukprev[n+1][prev]) # n+1 is the "current" state as we're going backwards, ukprev[n+1][prev] returns the state that maximized
        prev = ukprev[n+1][prev]
    backtrack = backtrack[::-1] # reverse the order
    return backtrack




'''
The forward algorithm calculates the probability of being in state i at time t, given all observations from time 1 to t. 

The forward algorithm is essentially identical to the Viterbi algorithm, with the major difference being that 
at each time step, we sum up all probabilities rather than choosing the maximizing transition.

The forward algorithm is also crucial for calculatating the likelihood of the observation given the entire model, i.e.

P(Observed Sequence | Model) = sum probabilities at the end point (last observation) over all states, yielding the total probability
of the observed sequence, given all current values of parameters. Sequences that better match the model will have higher probabilities
than those that do not.

Input: 
A - transition matrix
E - emission matrix
pi - initial state probabilities
states - list of states
alphabet - lets of letters comprising alphabet (e.g. ATCG)
O - observed sequence 

Output: Dictionary containing probabilities at each time point

'''

def fwd(O,A,pi,states,E):
    a = [{}]
    N = len(O)
    for state in states:
        a[0][state] = pi[state]+E[state][O[0]]
    for n in range(1,N):
        a.append({})
        if '$' in O[n]:
            for state in states:
                a[n][state] = a[n-1][state]
            continue
        for state in states:
            P=[]
            naiveP = []
            for pState in states:
                P.append(a[n-1][pState]+A[state][pState] + E[state][O[n]])
            P = logsumexp(P)
            a[n][state] = P
    return a
    
'''
Backward probabilties

Yields the probability of being in state i at time t, given all the observations for time t+1 to T,
where T is the last time point.

Similar to the foward algorithm, but moves in the opposite direction, considering the remainder of the sequence 
after the current observation.

Backward algorithm at time t = 1 will equal forward algorithm at time t = T

The forward and backward algorithms together calculate the 'smoothed probability' of being in state i at time t,
considering the entire sequence rather than prior and posterior information to the current observation. 


      Forward(t,state=i) * Backward(t,state=i)
---------------------------------------------------- = P(observation(t),state = i)
     SUM(Foward(t,state=i)*Backward(t,state=i))

Where the sum is taken over all states i in {query,null}
'''

def bkw(O,A,pi,states,E):
    b=[{}]
    N = len(O)

    for i in states:
        b[0][i] = 0

    for n in range(1,N):
        b.append({})
        if '$' in O[N-n]:
            for state in states:
                b[n][state] = b[n-1][state]
            continue
        for i in states:
            sumTerm = []
            for j in states:
                s = b[n-1][j]+A[i][j]+E[j][O[N-n]]
                sumTerm.append(s)
            P = logsumexp(sumTerm)
            b[n][i]= P
    b = b[::-1]
    return b

'''

Baum-Welch EM parameter updates

This is a custom implementation that only updates the transition matrix
the BW algorithm calculates the expected number of "observations" of each hidden state i, as well as the expected number of "observed"
transitions from each hidden state i to all possible states j

'''
def update(a,b,O,states,A,E):
    gamma = [{}]
    epsilon = [{'+':{'+':0,'-':0},'-':{'+':0,'-':0}}]
    T = len(O)
    for t in range(T-1):
        gamma.append({})
        epsilon.append({'+':{'+':0,'-':0},'-':{'+':0,'-':0}})
        # Fix for appended multiple training sequences
        if ('$' in O[t]) or ('$' in O[t+1]):
            for state in states:
                gamma[t][state] = gamma[t-1][state]
                for jstate in states:
                    epsilon[t][state][jstate] = epsilon[t-1][state][jstate] 
        else:
            norm = logsumexp([a[t]['+']+b[t]['+'],a[t]['-']+b[t]['-']])
            for i in states:
                gamma[t][i]=a[t][i]+b[t][i]-norm
                for j in states:
                    numerator = a[t][i]+A[i][j]+b[t+1][j]+E[j][O[t+1]]
                    denom = []
                    for k in states:
                        for w in states:
                            denom.append(a[t][k]+A[k][w]+b[t+1][w]+E[w][O[t+1]])
                    denom = logsumexp(denom)
                    epsilon[t][i][j] =numerator-denom
    margGamma = []
    margEpsilon = []
    for i in states:
        for j in states:
            for t in range(T-1):
                margGamma.append(gamma[t][i])
                margEpsilon.append(epsilon[t][i][j])
            A[i][j] = logsumexp(margEpsilon)-logsumexp(margGamma)
            margGamma = []
            margEpsilon = []
    return A
