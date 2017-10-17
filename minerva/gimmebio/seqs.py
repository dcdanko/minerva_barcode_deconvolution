


################################################################################

def baseToInt(base):
    if base == 'A':
        return 0
    elif base == 'C':
        return 1
    elif base == 'G':
        return 2
    elif base == 'T':
        return 3
    else:
        return 4

def rcBase(base):
    if base == 'A':
        return 'T'
    elif base =='C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'T':
        return 'A'
    else:
        return 'N'


def reverseComplement(kmer):
    rc = ''
    for base in kmer[::-1]:
        rc += rcBase( base)
    return rc
        
def canonical( kmer):
    rc = ''
    for base in kmer[::-1]:
        rc += rcBase( base)
        if kmer < rc:
            return kmer
    return rc
