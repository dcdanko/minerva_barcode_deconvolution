import sys
import argparse as ap
from .barcode_kmer_table import parseBarcodes, parseBarcodesAndRemoveStopKmers
from .build_and_filter_table import buildAndFilterTable
from .cluster_matrix import clusterDistMatrix
from .progress_bar import ProgressBar

################################################################################
#
# MAIN
#
################################################################################


def main():
    args = parseArgs()

    if args.remove_stop_kmers:
        barcodeTables = parseBarcodesAndRemoveStopKmers(sys.stdin,
                                                        args.K,
                                                        args.W,
                                                        args.dropout,
                                                        verbose=True)
    else:
        barcodeTables = parseBarcodes(sys.stdin, args.K, args.W, args.dropout,
                                      verbose=True)
    msg = '{:,} barcodes were at or above dropout threshold'
    print(msg.format(len(barcodeTables)), file=sys.stderr)

    totalAnchors = 0
    for bcTbl in barcodeTables:
        if bcTbl.numReads() >= args.anchor_dropout:
            totalAnchors += 1
    msg = '{:,} barcodes are anchors'
    print(msg.format(totalAnchors), file=sys.stderr)
    if totalAnchors == 0:
        print('No anchors found', file=sys.stderr)
        return
    progressBar = ProgressBar(totalAnchors)
    sys.stderr.write('\n')
    progressBar.write()

    output_file = sys.stdout
    if args.output_file != '-':
        output_file = open(args.output_file, 'w')
    anchors = (bcTbl for bcTbl in barcodeTables if bcTbl.numReads() >= args.anchor_dropout)
    for anchorTable in anchors:
        # build  and filter table
        sys.stderr.write(anchorTable.barcode)
        anchorTable = buildAndFilterTable(anchorTable, barcodeTables, args)
        if anchorTable is None:  # table was too small after filtering
            #            progressBar.increment()
            sys.stderr.write('\n')
            continue

        # reduce and cluster the rows
        readAssignments = clusterDistMatrix(anchorTable, args)
        writeClusters(anchorTable, readAssignments, output_file, args)
        sys.stderr.write('\n')
        #       progressBar.increment()
    if output_file != sys.stdout:
        output_file.close()


def writeClusters(anchorTable, readAssignments, output_file, args):
    for readId, clustNum in readAssignments.items():
        msg = '{}\t{}\t{}\n'.format(anchorTable.barcode, readId, clustNum)
        output_file.write(msg)


################################################################################
#
# ARGS
#
################################################################################


def parseArgs():
    parser = ap.ArgumentParser()

    parser.add_argument('-k', '--kmer-lens', dest='K',  default=20, type=int,
                        help='Lengths of kmers')
    parser.add_argument('-w', '--window-len', dest='W', default=40, type=int,
                        help='Window for sparse kmers')

    parser.add_argument('-d', '--dropout', dest='dropout', default=10, type=int,
                        help='Ignore barcodes with fewer reads')
    parser.add_argument('-a', '--anchor-dropout', dest='anchor_dropout', default=50, type=int,
                        help='Do not process anchors with fewer reads')

    parser.add_argument('--min-kmer', dest='rp_low_filter', default=1, type=float,
                        help='Filter kmers that rarely occur in other barcodes')
    parser.add_argument('--min-kmer-read', dest='min_kmer_per_read', default=1, type=float,
                        help='Require reads to have multiple kmers to overlap')
    parser.add_argument('--max-kmer', dest='rp_high_filter', default=0.03, type=float,
                        help='Filter kmers that occur constantly in other barcodes'
    parser.add_argument('--min-barcode', dest='bc_low_filter', default=2, type=float,
                        help='Filter barcodes with low kmer overlap')
    parser.add_argument('--max-barcode', dest='bc_high_filter', default=0.9, type=float,
                        help='Filter barcodes with high kmer overlap')

    parser.add_argument('--min-rows', dest='min_rows', default=5, type=int,
                        help='Do not process tables with fewer rows (kmers)')
    parser.add_argument('--min-cols', dest='min_cols', default=3, type=int,
                        help='Do not process tables with fewer cols (barcodes)')

    parser.add_argument('--eps',dest='dbscan_eps', default=0.26, type=float,
                        help='Distance threshold for DBSCAN clustering')
    parser.add_argument('--min-samples',dest='dbscan_min_samples',default=2, type=int,
                        help='Minimum samples in a cluster for dbscan')

    parser.add_argument('--output',dest='output_file',default='-', type=str,
                        help='File where results should be written to. Defaults to stdout.')


    # experimental args
    parser.add_argument('--remove-stopwords',dest='remove_stop_kmers', action='store_true',
                        help='Remove all kmers that occur 10x more often than average')


    parser.add_argument('--rescue-unassigned',dest='rescue_unassigned', action='store_true',
                        help='Try to assign reads not assigned by DBSCAN')

    parser.add_argument('--min-rescue', dest='min_rescue', default=2, type=int,
                        help='Rescue a read if it has this many connections to at most one cluster')

    parser.add_argument('--crack-edges', dest='crack_thresh', default=1, type=int,
                        help='Require that nodes that share a single edge share neighbours')


    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
