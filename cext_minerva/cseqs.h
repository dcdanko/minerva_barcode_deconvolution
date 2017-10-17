#include <Python.h>
#include <inttypes.h>


char complement( char base);

char * reverseComplement(char * seq, int seqLen);

char * canonicalInPlace(char * seq, int seqLen);

uint64_t nextRabinFingerprint( uint64_t oldHash, uint8_t oldBase, uint8_t newBase, uint8_t power, uint8_t radix);

PyObject * rabinFingerprints(const char * seq, uint8_t seqLen, uint8_t power, uint8_t radix);

PyObject * makeKmers(const char * seq, int seqLen, int k);
PyObject * makeCanonicalKmers(const char * seq, int seqLen, int k);
  
uint8_t baseToInt(char base);
