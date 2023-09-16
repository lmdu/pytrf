#define PY_SSIZE_T_CLEAN
#ifndef PYTRF_STR_H
#define PYTRF_STR_H
#include <Python.h>

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

} pytrf_STRFinder;

extern PyTypeObject pytrf_STRFinderType;

#endif
