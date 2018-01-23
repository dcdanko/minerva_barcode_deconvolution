import sys
import argparse as ap
import numpy as np
import math
from ..gimmebio.kmers import *


################################################################################
    
def passesFilter(nHits, possibleHits, minHits, maxHits):
    return belowFilter(nHits, possibleHits, maxHits) and aboveFilter(nHits, possibleHits, minHits)

def belowFilter(nHits, possibleHits, maxHits):
    if maxHits >= 1:
        return nHits <= maxHits
    else:
        pHits = nHits / possibleHits
        return pHits <= maxHits
    
def aboveFilter(nHits, possibleHits, minHits):
    if minHits >= 1:
        return nHits >= minHits
    else:
        pHits = nHits / possibleHits
        return pHits >= minHits

