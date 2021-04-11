#define PY_SSIZE_T_CLEAN
#ifndef STRIPY_VNTR_H
#define STRIPY_VNTR_H
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
	unsigned int min_motif;

	//max motif length
	unsigned int max_motif;

	//min repeats
	unsigned int min_repeats;

} stripy_VNTRMiner;

typedef struct {
	PyObject_HEAD

	//sequence name or identifier
	PyObject* seqid;

	//SSR motif sequence
	char* motif;

	//SSR start position
	Py_ssize_t start;

	//SSR stop position
	Py_ssize_t end;

	//number of tandem repeats
	unsigned int repeats;

	//SSR length
	unsigned int length;

} stripy_VNTR;

extern PyTypeObject stripy_VNTRMinerType;
extern PyTypeObject stripy_VNTRType;

#endif
