from itertools import product

# Cython code

# Read in data, type variables, create dictionary of kmers
# Set pseudo count to 1 for ALL kmers
# Call count function
# Return value of k and dictionary of k-mer counts

cpdef main(str fString,int k,str a):
    cdef unsigned int K = k
    cdef list kmers = [''.join(p) for p in product(a,repeat=K)]
#    cdef unsigned int pos = 4**K
    cdef unsigned int j
    cdef dict kmerDict = {}
    kmerDict = dict.fromkeys(kmers, 1)
#    for j in range(pos):
#        kmerDict[kmers[j]] = 1
    cdef unsigned long long N = len(fString)
    r = countKmers(fString,N,K,kmerDict)
    return k,r

# If k-mer is in dict of k-mers, add 1 to count
# Return k-mer count dict
cdef countKmers(str fString,unsigned long long N,unsigned int K,dict kmerDict):
    cdef unsigned long long i
    for i in range(0,N-K+1):
        curr = fString[i:i+K]
        if curr in kmerDict:
            kmerDict[curr]+=1
        else:
          continue
    return kmerDict
