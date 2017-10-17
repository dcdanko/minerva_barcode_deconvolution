from sklearn.cluster import dbscan

################################################################################
#
# Reductions and clustering
#
################################################################################

def clusterDistMatrix( anchorTable, args):
    kmerAssignments = {}    
    iam, rowNames = anchorTable.asInverseAdjacency()

    _, labels = dbscan( iam,
                        eps=args.dbscan_eps,
                        min_samples=args.dbscan_min_samples,
                        metric='precomputed')
    
    for ind, kmer in rowNames.items():
        clust = labels[ind]
        if clust == -1:
            continue
        kmerAssignments[kmer] = clust

    return kmerAssignments


