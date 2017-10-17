#include "cseqs.h"
#include <Python.h>

static char  module_docstring[] =
  "A module for sequence related functions";
static char reverseComplement_docstring[] =
  "Get the reverse complement of a sequence";
static char canonical_docstring[] =
  "Get the canonical version of a sequence";
static char nextRabinFingerprint_docstring[] =
  "Calculate the next rabin fingerprint";
static char rabinFingerprints_docstring[] =
  "Calculate rabin fingerprints for a sequence";
static char makeKmers_docstring[] =
  "Make kmers";
static char makeCanonicalKmers_docstring[] =
  "Make canonical kmers";
  

static PyObject * cseqs_nextRabinFingerprint(PyObject* self, PyObject* args);
static PyObject * cseqs_rabinFingerprints(PyObject* self, PyObject* args);
static PyObject * cseqs_makeKmers(PyObject *self, PyObject *args);
static PyObject * cseqs_makeCanonicalKmers(PyObject *self, PyObject *args);
static PyObject * cseqs_reverseComplement(PyObject* self, PyObject* args);
static PyObject * cseqs_canonical(PyObject* self, PyObject* args);

static PyMethodDef cseqs_methods[] = {
  {"makeKmers", cseqs_makeKmers, METH_VARARGS, makeKmers_docstring},
  {"makeCanonicalKmers", cseqs_makeCanonicalKmers, METH_VARARGS, makeCanonicalKmers_docstring},    
  {"nextRabinFingerprint", cseqs_nextRabinFingerprint, METH_VARARGS, nextRabinFingerprint_docstring},
  {"rabinFingerprints", cseqs_rabinFingerprints, METH_VARARGS, rabinFingerprints_docstring},
  {"reverseComplement", cseqs_reverseComplement, METH_VARARGS, reverseComplement_docstring},
  {"canonical", cseqs_canonical, METH_VARARGS, canonical_docstring}
};

static struct PyModuleDef cseqsmodule = {
  PyModuleDef_HEAD_INIT,
  "cseqs",
  module_docstring,
  -1,
  cseqs_methods
};

PyMODINIT_FUNC PyInit_cseqs(void){
  return PyModule_Create(&cseqsmodule);
}

static PyObject * cseqs_reverseComplement(PyObject *self, PyObject *args){
  char * seq;
  int seqLen;
  if (!PyArg_ParseTuple(args, "si", &seq, &seqLen)){
    return NULL;
  }
  char * rc = reverseComplement(seq, seqLen);
  PyObject *ret= Py_BuildValue("s", rc);
  return ret;
    
}

static PyObject * cseqs_canonical(PyObject *self, PyObject *args){
  char * seq;
  int seqLen;
  if (!PyArg_ParseTuple(args, "si", &seq, &seqLen)){
    return NULL;
  }
  char * canon = canonicalInPlace(seq, seqLen);
  PyObject * ret= Py_BuildValue("s", canon);
  return ret;
    
}

static PyObject * cseqs_makeKmers(PyObject *self, PyObject *args){
  const char * seqObj;
  int seqLen, k;

  if (!PyArg_ParseTuple(args, "sii", &seqObj, &seqLen, &k)){
    return NULL;
  }
  PyObject * ret = makeKmers(seqObj, seqLen, k);
  return ret;    
}

static PyObject * cseqs_makeCanonicalKmers(PyObject *self, PyObject *args){
  const char * seqObj;
  int seqLen, k;

  if (!PyArg_ParseTuple(args, "sii", &seqObj, &seqLen, &k)){
    return NULL;
  }
  PyObject * ret = makeCanonicalKmers(seqObj, seqLen, k);
  return ret;    
}


static PyObject * cseqs_nextRabinFingerprint(PyObject *self, PyObject *args){
  uint64_t oldHash;
  uint8_t oldBase, newBase, power, radix;
  if (!PyArg_ParseTuple(args, "Lbbbb", &oldHash, &oldBase, &newBase, &power, &radix)){
    return NULL;
  }

  uint64_t newHash = nextRabinFingerprint(oldHash, oldBase, newBase, power, radix);
  PyObject *ret= Py_BuildValue("L", newHash);
  return ret;
    
}

static PyObject * cseqs_rabinFingerprints(PyObject *self, PyObject *args){
  const char * seqObj;
  uint8_t seqLen, power, radix;
  if (!PyArg_ParseTuple(args, "sbbb", &seqObj, &seqLen, &power, &radix)){
    return NULL;
  }

  return rabinFingerprints(seqObj, seqLen, power, radix);
    
}


