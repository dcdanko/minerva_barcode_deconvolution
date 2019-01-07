import sys
import click
from .kmers import MinSparseKmerSet

################################################################################

class BarcodeKmerTable:

    def __init__(self, barcode):
        self.barcode = barcode
        self.nreads = 0
        self.readKmerSets = {}
        self.reverseMap = None
        self.kmersWithCounts = {}

    def addKmersFromRead(self, readId, kmerSet):
        self.nreads += 1
        self.readKmerSets[readId] = kmerSet
        for kmer in kmerSet:
            try:
                self.kmersWithCounts[kmer] += 1
            except KeyError:
                self.kmersWithCounts[kmer] = 1
                    
    def _buildReverseMap(self):
        self.reverseMap = {}
        for rp, kmerSet in self.readKmerSets.items():
            for kmer in kmerSet:
                try:
                    self.reverseMap[kmer].append(rp)
                except KeyError:
                    self.reverseMap[kmer] = [rp]

    def kmerToReadPairs(self, kmer):
        if self.reverseMap is None:
            self._buildReverseMap()
        return self.reverseMap[kmer]

    def allKmers(self):
        # it is important that kmers that occur multiple times are returned multiple times
        for kmer, count in self.kmersWithCounts.items():
            for _ in range(count):
                yield kmer

    def containsKmer(self, kmer):
        return kmer in self.kmersWithCounts

    def __iter__(self):
        return iter(self.allKmers())

    def __len__(self):
        return self.nreads
    
################################################################################
    
def iterBarcodesAsDocuments(filelike, k, w, dropout, verbose=False):
    for nBC, bc in enumerate( iterReadCloud(filelike)):
        if verbose and ((nBC % 100) == 0):
            sys.stderr.write('\rparsed {:,} barcodes'.format(nBC))
        if len(bc) < args.dropout:
            continue
        
        bcKmers = BarcodeKmerTable(bc.barcode)
        for rp in bc:
            kmerSet = MinSparseKmerSet( k, w, [rp.r1.seq, rp.r2.seq])
            bcKmers.addKmersFromRead(rp.sid, kmerSet)
        yield bcKmers
            
    if verbose:
        sys.stderr.write('\rparsed {:,} barcodes\n'.format(nBC))


@click.command()
@click.option('-k', default=20, help='Lengths of kmers')
@click.option('-w', default=40, help='Length of window (in bases, not kmers)')
@click.option('-d', '--dropout', default=100, help='Ignore smaller read clouds')
@click.option('-v/-q', '--verbose/--quiet', default=True, help='Periodically report output')
def printBarcodesAsDocuments(k, w, dropout, verbose):
    for bc in iterBarcodesAsDocuments(sys.stdin, k, w, dropout, verbose=verbose):
        sys.stdout.write(bc.barcode)
        for kmer in bc:
            sys.stdout.write('\t{}'.format(kmer))
        sys.stdout.write('\n')

################################################################################

def findStopKmers(filelike, k, w, dropout, meanMultiplier=10, verbose=False):
    kmerCounts = {}
    totalCount = 0
    for nBC, bc in enumerate( iterReadCloud(filelike)):
        if verbose and ((nBC % 100) == 0):
            sys.stderr.write('\rparsed {:,} barcodes'.format(nBC))
        if len(bc) < args.dropout:
            continue
        for rp in bc:
            kmerSet = MinSparseKmerSet( k, w, [rp.r1.seq, rp.r2.seq])
            for kmer in kmerSet:
                totalCount += 1
                try:
                    kmerCounts[kmer] += 1
                except KeyError:
                    kmerCounts[kmer] = 1
    if verbose:
        sys.stderr.write('\rparsed {:,} barcodes\n'.format(nBC))

    aveCount = totalCount / len(kmerCounts)
    cutoff = aveCount * meanMultiplier
    for kmer, count in kmerCoutns.items():
        if count > cutoff:
            yield kmer

            
        
@click.command()
@click.option('-k', default=20, help='Lengths of kmers')
@click.option('-w', default=40, help='Length of window (in bases, not kmers)')
@click.option('-m', default=10, help='Find kmers that are more than Nx as common as average')
@click.option('-d', '--dropout', default=100, help='Ignore smaller read clouds')
@click.option('-v/-q', '--verbose/--quiet', default=True, help='Periodically report output')
def printStopKmers(k, w, m, dropout, verbose):
    for kmer in findStopKmers(sys.stdin, k, w, dropout,
                            meanMultiplier=m,
                            verbose=verbose):
        print(kmer)

            
if __name__ == '__main__':
    printBarcodesAsDocuments()
