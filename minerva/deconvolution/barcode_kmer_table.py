from minerva.gimmebio.kmers import MinSparseKmerSet
from minerva.gimmebio.readclouds import iterReadClouds
from itertools import combinations
import numpy as np
import sys

################################################################################

def parseBarcodesAndRemoveStopKmers(filelike, K, W, dropout,
                                    meanMultiplier=10,
                                    verbose=False):
    bcTbls = parseBarcodes(filelike, K, W, dropout, verbose=verbose)
    kmerCounts = {}
    totalCount = 0
    for bcTbl in bcTbls:
        for kmer in bcTbl.allKmers():
            totalCount += 1
            try:
                kmerCounts[kmer] += 1
            except KeyError:
                kmerCounts[kmer] = 1
    aveCount = totalCount / len(kmerCounts)
    cutoff = aveCount * meanMultiplier
    stopKmers = set()
    for kmer, count in kmerCounts.items():
        if (count == 1) or (count >  cutoff):
            stopKmers.add(kmer)
    print('Removing {:,} stop and singleton kmers'.format( len(stopKmers)), file=sys.stderr)
    for bcTbl in bcTbls:
        bcTbl.removeKmers( stopKmers)
    print('Removed stop and singleton kmers', file=sys.stderr)
    return bcTbls

def parseBarcodes(filelike, K, W, dropout, verbose=False):
    bcTbls = []
    for nBC, bc in enumerate( iterReadClouds(filelike)):
        if verbose and ((nBC % 100) == 0):
            sys.stderr.write('\rparsed {:,} barcodes'.format(nBC))
        if len(bc) < dropout: 
            continue
        readKmerSets = {}
        for rp in bc:
            seqs = [rp.r1.seq, rp.r2.seq]
            readKmerSets[rp.sid] = MinSparseKmerSet( K, W, seqs)
        bcTbls.append( BarcodeTable(bc.barcode, readKmerSets))
    if verbose:
        sys.stderr.write('\rparsed {:,} barcodes\n'.format(nBC))
    return bcTbls

        
################################################################################


class BarcodeTable:

    def __init__(self, barcode, readKmerSets):
        self.barcode = barcode
        self.readKmerSets = readKmerSets
        self.reverseMap = None
        self.tbl = {}
        self._allKmers = set()
        for rp, kmerSet in self.readKmerSets.items():
            for kmer in kmerSet:
                self._allKmers.add(kmer)

    def _buildReverseMap(self):
        self.reverseMap = {}
        for rp, kmerSet in self.readKmerSets.items():
            for kmer in kmerSet:
                try:
                    self.reverseMap[kmer].append(rp)
                except KeyError:
                    self.reverseMap[kmer] = [rp]

    def removeKmers(self, stopKmers):
        assert self.reverseMap is None
        self._allKmers = {k for k in self._allKmers if k not in stopKmers}
        for rp, kmerSet in self.readKmerSets.items():
            kmerSet.removeKmers( stopKmers)
        
                   
    ############################################################################


    def hasColumn(self, otherTbl):
        return (otherTbl.barcode in self.tbl) or (otherTbl.barcode == self.barcode)

    def kmerToReadPairs(self, kmer):
        if not self.reverseMap:
            self._buildReverseMap()
        return self.reverseMap[kmer]

    def setColumn(self, barcode, readPairs):
        self.tbl[barcode] = readPairs

    def kmerSet(self):
        return self._allKmers
        
    def allKmers(self):
        # it is important that kmers that occur multiple times are returned multiple times
        for kmerSet in self.readKmerSets.values():
            for kmer in kmerSet:
                yield kmer

    def containsKmer(self, kmer):
        return kmer in self._allKmers


    def asAdjacencyList(self):
        rpInds = {}
        rowNames = {}
        for bc, rps in self.tbl.items():
            for rp in rps:
                if rp not in rpInds:
                    n = len(rpInds)
                    rpInds[rp] = n
                    rowNames[n] = rp
                    
        adjList = {}
        for bc, rps in self.tbl.items():
            for r1, r2 in combinations(rps, 2):
                try:
                    adjList[r1].add(r2)
                except KeyError:
                    adjList[r1] = set([r2])

                try:
                    adjList[r2].add(r1)
                except KeyError:
                    adjList[r2] = set([r1])
        return adjList

    
    def asInverseAdjacency(self, crackThresh=-1):
        rpInds = {}
        rowNames = {}
        for bc, rps in self.tbl.items():
            for rp in rps:
                if rp not in rpInds:
                    n = len(rpInds)
                    rpInds[rp] = n
                    rowNames[n] = rp
                    
        ws = {}
        for bc, rps in self.tbl.items():
            for r1, r2 in combinations(rps, 2):
                key = tuple( sorted( [r1,r2]))
                try:
                    ws[key] += 1
                except KeyError:
                    ws[key] = 1

        if crackThresh >= 1:
            ws = self._crackEdges(ws, crackThresh)
                    
        iam = (1000*1000) * np.ones( (len(rpInds), len(rpInds)))
        for (r1,r2), w in ws.items():
            i1 = rpInds[r1]
            i2 = rpInds[r2]
            iam[i1,i2] = 1.0 / w
            iam[i2,i1] = 1.0 / w            

        sys.stderr.write(' -> INV_ADJ:' + str(iam.shape))                    
        return iam, rowNames


    def _crackEdges(self, ws, crackThresh):
        out = {}
        adjList = self.asAdjacencyList()
        for (r1,r2), w in ws.items():
            if w == 1:
                sharedNeibs = adjList[r1] & adjList[r2]
                if len(sharedNeibs) < crackThresh:
                    continue
            out[(r1,r2)] = w
        return out
    
    def asDataFrame(self):
        df = {}
        for bc, rps in self.tbl.items():
            df[bc] = {rp:1 for rp in rps}
        df = pd.DataFrame(df)
        df.fillna(value=0, inplace=True)
        return df

    def asDict(self):
        return self.tbl

    def numColumns(self):
        return len(self.tbl)

    def numReads(self):
        return len(self.readKmerSets)

    def numRows(self):
        uniqueRows = set()
        for rps in self.tbl.values():
            for rp in rps:
                uniqueRows.add(rp)
        return len(uniqueRows)

    def shape(self):
        return self.numRows(), self.numColumns()
