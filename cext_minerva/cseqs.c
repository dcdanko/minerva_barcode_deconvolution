#include <Python.h>
#include "cseqs.h"
#include <math.h>
#include <stdio.h>



char complement(char base){
  switch(base){
  case 'A':
    return 'T';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'T':
    return 'A';
  default:
    return 'N';
  }
}

char * reverseComplement(char * seq,  int seqLen){
  char * rc = malloc((seqLen+1)*sizeof(char));
  int i;
  for(i=0; i<seqLen; i++){
    rc[i] = complement( seq[seqLen - 1 - i]);
  }
  rc[seqLen] = '\0';
  return rc;
}


char * canonicalInPlace(char * seq, int seqLen){
  char * rc = reverseComplement( seq, seqLen);
  int i, j;
  for(i=0; i<seqLen; i++){
    if(seq[i] < rc[i]){
      free(rc);
      return seq;
    } else if(seq[i] > rc[i]){
      for(j=0; j<seqLen; j++){
	seq[j] = rc[j];
      }
      free(rc);
      return seq;
    }
  }
  free(rc);
  return seq;
}

uint64_t nextRabinFingerprint( uint64_t oldHash, uint8_t oldBase, uint8_t newBase, uint8_t power, uint8_t radix){
  return radix * (oldHash - (oldBase * (radix << power))) + newBase;
}

uint8_t baseToInt(char base){
  switch(base){
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return 4;
  }
}

PyObject * makeKmers(const char * seq, int seqLen, int k){
  int nkmers;
  nkmers = seqLen - k + 1;
  int i,j;

  PyObject * ret = PyList_New(nkmers);
  char kmer[k+1];
  for(i=0; i<nkmers; i++){
    for(j=0; j<k; j++){
      kmer[j] = seq[i+j];
    }
    kmer[k] = '\0';
    PyList_SET_ITEM(ret, i, Py_BuildValue("s", kmer));    
  }
  return ret;
}

PyObject * makeCanonicalKmers(const char * seq, int seqLen, int k){
  int nkmers;
  nkmers = seqLen - k + 1;
  int i,j;

  PyObject * ret = PyList_New(nkmers);
  char kmer[k+1];
  for(i=0; i<nkmers; i++){
    for(j=0; j<k; j++){
      kmer[j] = seq[i+j];
    }
    kmer[k] = '\0';
    PyList_SET_ITEM(ret, i, Py_BuildValue("s", canonicalInPlace(kmer, k)));    
  }
  return ret;
}


PyObject * rabinFingerprints(const char * seq, uint8_t seqLen, uint8_t power, uint8_t radix){
  uint8_t useq[seqLen];
  int i, k;
  uint64_t initHash;
  for(i=0; i<seqLen; i++){
    useq[i] = baseToInt( seq[i]);
  }
  
  k = 1 + pow(2, power);
  initHash = 0;
  uint64_t allHash[seqLen - k + 1];
  for(i=0; i<k; i++){
    initHash += useq[i] * pow(radix, (k - 1 - i));
  }
  allHash[0] = initHash;

  for(i=k; i<seqLen; i++){
    allHash[i - k + 1] = nextRabinFingerprint( allHash[i - k],
					       useq[i - k],
					       useq[i],
					       power,
					       radix);

  }

  char sub[k+1];
  int j;
  PyObject * ret = PyList_New(seqLen - k + 1);
  for(i=0; i<(seqLen - k + 1); i++){
    PyObject * tup = PyTuple_New(2);
    for(j=i; j<i+k; j++){
      sub[j-i] = seq[j];
    }
    sub[k] = '\0';
    PyTuple_SET_ITEM( tup, 0, Py_BuildValue("s", sub));
    PyTuple_SET_ITEM( tup, 1, PyLong_FromLongLong( allHash[i]));    
    PyList_SET_ITEM(ret, i, tup);
  }
  return ret;
}
  
  
