


################################################################################

def iterBarcodes(filelike, K, W, dropout, verbose=False):
    for nBC, bc in enumerate( iterReadCloud(filelike)):
        if verbose and ((nBC % 100) == 0):
            sys.stderr.write('\rparsed {:,} barcodes'.format(nBC))
        if len(bc) < dropout: 
            continue
        readKmerSets = {}
        for rp in bc:
            seqs = [rp.r1.seq, rp.r2.seq]
            readKmerSets[rp.sid] = MinSparseKmerSet( K, W, seqs)
        yield BarcodeTable(bc.barcode, readKmerSets)
    if verbose:
        sys.stderr.write('\r{:,} read clouds\n'.format(nBC))
    
################################################################################


class BarcodeTable:

    def __init__(self, barcode, readKmerSets, isAnchor):
        self.barcode = barcode
        self.numReads = len(readKmerSets)
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

    def hasColumn(self, otherTbl):
        return (otherTbl.barcode in self.tbl) or (otherTbl.barcode == self.barcode)

    def kmerToReadPairs(self, kmer):
        if not self.reverseMap:
            self._buildReverseMap()
        return self.reverseMap[kmer]

    def setColumn(self, barcode, readPairs):
        self.tbl[barcode] = readPairs

    def allKmers(self):
        if self.
        # it is important that kmers that occur multiple times are returned multiple times
        for kmerSet in self.readKmerSets.values():
            for kmer in kmerSet:
                yield kmer

    def containsKmer(self, kmer):
        return kmer in self._allKmers


    def asInverseAdjacency(self):
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

        iam = (1000*1000) * np.ones( (len(rpInds), len(rpInds)))
        for (r1,r2), w in ws.items():
            i1 = rpInds[r1]
            i2 = rpInds[r2]
            iam[i1,i2] = 1.0 / w
            iam[i2,i1] = 1.0 / w            

        sys.stderr.write(' -> INV_ADJ:' + str(iam.shape))                    
        return iam, rowNames

    
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
