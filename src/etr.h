/* 
 * etr.h -- Exact Tandem Repeat
 *
 * define the base class for exact tandem repeat element
 *
 */

#ifndef PYTRF_ETR_H
#define PYTRF_ETR_H
#include "Python.h"

typedef struct {
	PyObject_HEAD

	//sequence name or identifier
	PyObject* seqid;

	//tandem repeat element motif sequence
	char* motif;

	//tandem repeat element motif length
	int mlen;

	//tandem repeat element start position
	Py_ssize_t start;

	//tandem repeat element stop position
	Py_ssize_t end;

	//number of tandem repeats
	int repeats;

	//tandem repeat element length
	int length;

} pytrf_ETR;

extern PyTypeObject pytrf_ETRType;

#endif
