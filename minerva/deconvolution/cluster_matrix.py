from sklearn.cluster import dbscan

################################################################################
#
# Reductions and clustering
#
################################################################################

def clusterDistMatrix( anchorTable, args):
    kmerAssignments = {}    
    iam, rowNames = anchorTable.asInverseAdjacency(crackThresh=args.crack_thresh)

    _, labels = dbscan( iam,
                        eps=args.dbscan_eps,
                        min_samples=args.dbscan_min_samples,
                        metric='precomputed')
    
    for ind, kmer in rowNames.items():
        clust = labels[ind]
        if clust == -1:
            continue
        kmerAssignments[kmer] = clust

    if args.rescue_unassigned:
        rescueUnassigned( kmerAssignments, anchorTable, args)
        
    return kmerAssignments


def rescueUnassigned(assignments, anchorTable, args):
    adjacency = anchorTable.asAdjacencyList()
    unassigned = [kmer for kmer in adjacency.keys() if kmer not in assignments]

    clusters = {}
    for kmer in unassigned:
        for neib in adjacency[kmer]:
            try:
                clust = assignments[neib]
            except KeyError:
                continue

            try:
                clusters[kmer][clust] += 1
            except KeyError:
                try:
                    clusters[kmer][clust] = 1
                except KeyError:
                    clusters[kmer] = { clust : 1}

    for kmer, cls in clusters.items():
        potentials = []
        for cl, count in cls.items():
            if count >= args.min_rescue:
                potentials.append(cl)
        if len(potentials) == 1:
            assignments[kmer] = potentials[0]


            
            
    
