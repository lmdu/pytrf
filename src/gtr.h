/* 
 * gtr.h -- Generic tandem repeat header
 *
 * define the finder class for generic tandem repeat
 *
 */

#define PY_SSIZE_T_CLEAN
#ifndef PYTRF_GTR_H
#define PYTRF_GTR_H
#include <Python.h>

typedef struct {
	PyObject_HEAD

	//input sequence name
	PyObject *seqname;

	//input sequence object
	PyObject *seqobj;

	//sequence length
	Py_ssize_t size;

	//next start position for tandem repeat identification
	Py_ssize_t next_start;

	//pointer to sequence object
	const char *seq;

	//max motif length
	int max_motif;

	//min repeats
	int min_repeat;

	//min length
	int min_length;

} pytrf_GTRFinder;

extern PyTypeObject pytrf_GTRFinderType;

#endif
