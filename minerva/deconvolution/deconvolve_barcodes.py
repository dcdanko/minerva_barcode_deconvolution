import sys
import argparse as ap
from .barcode_kmer_table import iterBarcodes
from .build_and_filter_table import buildAndFilterTable
from .cluster_matrix import clusterDistMatrix

################################################################################
#
# MAIN
#
################################################################################

def main():
    args = parseArgs()
    
    barcodeIter = iterBarcodes( sys.stdin, args.K, args.W, args.dropout, verbose=True)
    barcodeTables = [bcTbl for bcTbl in barcodeIter]

    msg = '{:,} barcodes were at or above dropout threshold'
    print(msg.format(len(barcodeTables)), file=sys.stderr)

    anchors = (bcTbl for bcTbl in barcodeTables if bcTbl.numReads > args.anchor_dropout)
    for  anchorTable in anchors:
        sys.stderr.write('\n'+anchorTable.barcode)

        # build  and filter table
        anchorTable = buildAndFilterTable(anchorTable, barcodeTables, args)
        if anchorTable is None: # table was too small after filtering
            continue

        # reduce and cluster the rows
        readAssignments = clusterDistMatrix( anchorTable, args)

        writeClusters( anchorTable, readAssignments, args)


        
def writeClusters(anchorTable, readAssignments, args):
    for readId, clustNum in readAssignments.items():
        msg = '{}\t{}\t{}'.format(anchorTable.barcode, readId, clustNum)
        print(msg)


################################################################################
#
# ARGS
#
################################################################################
        
def  parseArgs():
    parser = ap.ArgumentParser()

    parser.add_argument('-k', '--kmer-lens', dest='K',  default=20, type=int,
                        help='Lengths of kmers')
    parser.add_argument('-w', '--window-len', dest='W', default=40, type=int,
                        help='Window for sparse kmers')
    
    parser.add_argument('-d', '--dropout', dest='dropout', default=100, type=int,
                        help='Ignore barcodes with fewer reads')
    parser.add_argument('--anchor-dropout', dest='anchor_dropout', default=200, type=int,
                        help='Do not process anchors with fewer reads')

    parser.add_argument('--min-kmer', dest='rp_low_filter', default=1, type=float,
                        help='Filter kmers that rarely occur in other barcodes')
    parser.add_argument('--min-kmer-read', dest='min_kmer_per_read', default=1, type=float,
                        help='Require reads to have multiple kmers to overlap')
    parser.add_argument( '--max-kmer', dest='rp_high_filter', default=0.5, type=float,
                         help='Filter kmers that occur constantly in other barcodes')    
    parser.add_argument('--min-barcode', dest='bc_low_filter', default=2, type=float,
                        help='Filter barcodes with low kmer overlap')
    parser.add_argument('--max-barcode', dest='bc_high_filter', default=0.9, type=float,
                        help='Filter barcodes with high kmer overlap')    
    
    parser.add_argument('--min-rows', dest='min_rows', default=5, type=int,
                        help='Do not process tables with fewer rows (kmers)')    
    parser.add_argument('--min-cols', dest='min_cols', default=3, type=int,
                        help='Do not process tables with fewer cols (barcodes)')

    parser.add_argument('--dbscan-eps',dest='dbscan_eps', default=0.3, type=float,
                        help='Distance threshold for DBSCAN clustering')
    parser.add_argument('--dbscan-min-samples',dest='dbscan_min_samples',default=3, type=int,
                        help='Minimum samples in a cluster for dbscan')    

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
