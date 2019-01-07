
from .seqs import *
import click

################################################################################

class MinSparseKmerSet:
    '''
    Represents the kmers in a read
    '''
    def __init__(self, K, W, seqs, canonical=True):
        self.K = K
        self.W = W
        self.nseqs = len(seqs)
        self.makeKmers(seqs, canonical)
        
    def makeKmers(self, seqs, canonical):
        self.kmers = {}
        for seq in seqs:
            kmers = makeKmers(seq, self.K, canon=canonical)
            kmers = [(kmer, hash(kmer)) for kmer in kmers]
            kmersInWindow = self.W - self.K + 1
            numWindows =  len(kmers) - kmersInWindow + 1
            for windowStart in range( numWindows):
                window = kmers[windowStart:windowStart + kmersInWindow]

                if windowStart == 0:
                    # first window, minimum does not exist yet
                    minKmer = min(window,key=lambda x: x[1])
                    self.kmers[minKmer[0]] = 1 
                elif kmers[windowStart - 1][1] == minKmer[1]:
                    # the minimum has dropped out of the window
                    minKmer = min(window,key=lambda x: x[1])
                    self.kmers[minKmer[0]] = 1                     
                elif window[-1][1] < minKmer[1]:
                    # the new kmer is lower
                    minKmer = window[-1]
                    self.kmers[minKmer[0]] = 1                     
                    
                
                
                
    def getCount(self, kmer):
        try:
            return self.kmers[kmer]
        except KeyError:
            return 0

    def overlap(self, other):
        'Return a list of kmers that occur in both sets'
        
        assert type(other) == type(self)
        out = []
        for kmer in other:
            if kmer in self:
                out.append(kmer)
        return out

    def remove(self, kmer):
        del self.kmers[kmer]

    def removeKmers(self, stopKmers):
        self.kmers = {k:v for k,v in self.kmers.items() if k not in stopKmers}
        
    def __contains__(self, kmer):
        return kmer in self.kmers

    def __iter__(self):
        return iter(self.kmers.keys())

    def withCounts(self):
        return self.kmers.items()
    
    def __str__(self):
        out = sorted( self.kmers.keys())
        return ' '.join(out)

    def __len__(self):
        return len(self.kmers)




def makeKmers(seq, K, canon=True):
    out = []
    for start in range(len(seq)-K+1):
        kmer = seq[start:start+K]
        if canon:
            kmer = canonical(kmer)
        out.append(kmer)
    
    return out
