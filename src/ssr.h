#define PY_SSIZE_T_CLEAN
#ifndef STRIPY_SSR_H
#define STRIPY_SSR_H
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

	//minimal length for mono-, di-, tri-, tetra-, penta- and hexa- types
	Py_ssize_t min_lens[7];

} stripy_SSRMiner;

typedef struct {
	PyObject_HEAD

	//sequence name or identifier
	PyObject* seqid;

	//SSR motif sequence
	char motif[7];

	//SSR start position
	Py_ssize_t start;

	//SSR stop position
	Py_ssize_t end;

	//number of tandem repeats
	unsigned int repeats;

	//SSR length
	unsigned int length;

} stripy_SSR;

extern PyTypeObject stripy_SSRMinerType;
extern PyTypeObject stripy_SSRType;

#endif
