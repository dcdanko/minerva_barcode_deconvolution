from .utils import passesFilter
from collections import Counter
import sys


################################################################################
#
# Build Table
#
################################################################################

def buildAndFilterTable( anchorTable, barcodeTables, args):
    # build and filter table
    buildTable(anchorTable, barcodeTables, args)
    if anchorTable.numColumns() == 0:
        sys.stderr.write('\tempty_table_before_filter')
        return None
    filterTable( anchorTable, args)

    nrow, ncol = anchorTable.shape()
    tableTooSmall = nrow < args.min_rows or ncol < args.min_cols
    tableIsEmpty = nrow ==0 or ncol == 0

    if not tableIsEmpty:
        sys.stderr.write(' -> RF:' + str(anchorTable.shape()))
    else:
        sys.stderr.write('\t empty_table')

    if tableTooSmall or tableIsEmpty:
        return None

    return anchorTable


def buildTable(anchorTable, otherTables, args):
    for otherTable in otherTables:
        olap = len(otherTable.kmerSet() & anchorTable.kmerSet())
        maxOlap = anchorTable.numReads()
        if passesFilter( olap,
                         maxOlap,
                         args.bc_low_filter,
                         args.bc_high_filter):
            
            if not anchorTable.hasColumn( otherTable):
                buildNewColumns(anchorTable, otherTable, args)
            
def buildNewColumns(anchorTable, otherTable, args):
    kmerCounts = {}
    for anchorKmer in anchorTable.allKmers():
        if otherTable.containsKmer(anchorKmer):
            try:
                kmerCounts[anchorKmer] += 1
            except KeyError:
                kmerCounts[anchorKmer] = 1
    validKmers = [kmer for kmer, count in kmerCounts.items() if count == 1]
    anchorReadPairs = {}
    for kmer in validKmers:
        rp = anchorTable.kmerToReadPairs(kmer)[0] 
        try:
            anchorReadPairs[rp] += 1
        except KeyError:
            anchorReadPairs[rp] = 1
    anchorReadPairs = [rp for rp, count in anchorReadPairs.items() if count >= args.min_kmer_per_read]
    anchorTable.setColumn(otherTable.barcode, anchorReadPairs)
    if otherTable.numReads() >= args.anchor_dropout:
        otherReadPairs = {}
        for kmer in validKmers:
            rp = otherTable.kmerToReadPairs(kmer)[0] 
            try:
                otherReadPairs[rp] += 1
            except KeyError:
                otherReadPairs[rp] = 1
        otherReadPairs = [rp for rp, count in otherReadPairs.items() if count >= args.min_kmer_per_read]

        otherTable.setColumn(anchorTable.barcode, otherReadPairs)
        
def filterTable(anchorTbl, args):
    sys.stderr.write(' RU:({}, {})'.format(anchorTbl.numReads(), anchorTbl.numColumns()))    
    
    filtBarcodeTbl = {}
    for barcode, readPairs in anchorTbl.asDict().items():
        numReadPairs = len(readPairs)
        if passesFilter( numReadPairs,
                         anchorTbl.numReads(),
                         args.bc_low_filter,
                         args.bc_high_filter):
            filtBarcodeTbl[barcode] = readPairs

    rpCounts = Counter()
    nbc  = 0
    for barcode, readPairs in filtBarcodeTbl.items():
        nbc += 1
        rpCounts.update(readPairs)

    sys.stderr.write(' -> bcFilt:({}, {})'.format(len(rpCounts), nbc))    
        
    goodReadPairs = set()
    for rpSID, numBarcodes in rpCounts.items():
        if passesFilter( numBarcodes,
                         len(filtBarcodeTbl),
                         args.rp_low_filter,
                         args.rp_high_filter):
            goodReadPairs.add(rpSID)
        
    filtTable = {}
    for barcode, readPairs in filtBarcodeTbl.items():
        filtRps = []
        for rpSID in readPairs:
            if rpSID in goodReadPairs:
                filtRps.append(rpSID)
        filtTable[barcode] = filtRps

    anchorTbl.tbl = filtTable

