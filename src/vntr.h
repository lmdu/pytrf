#define PY_SSIZE_T_CLEAN
#ifndef STRIA_VNTR_H
#define STRIA_VNTR_H
#include "Python.h"

typedef struct {
	PyObject_HEAD

	//input sequence name
	PyObject *seqname;

	//input sequence object
	PyObject *seqobj;

	//pointer to sequence object
	const char *seq;

	//sequence length
	Py_ssize_t size;

	//next start position for tandem repeat identification
	Py_ssize_t next_start;

	//minimal motif length
	int min_motif;

	//max motif length
	int max_motif;

	//min repeats
	int min_repeat;

} stria_VNTRMiner;

extern PyTypeObject stria_VNTRMinerType;

#endif
