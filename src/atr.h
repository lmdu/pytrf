/*
 * atr.h -- approximate or imperfect tandem repeat element
 *
 */

#define PY_SSIZE_T_CLEAN
#ifndef PYTRF_ATR_H
#define PYTRF_ATR_H
#include "Python.h"

typedef struct {
	PyObject_HEAD

	//input sequence name
	PyObject *seqid;

	//input sequence obj
	PyObject *seqobj;

	//motif sequence
	PyObject *motif;

	//imperfect tandem repeat start position
	Py_ssize_t start;

	//imperfect tandem repeat stop position
	Py_ssize_t end;

	//motif length
	int mlen;

	//tandem length
	int length;

	//number of matches
	int matches;

	//number of substitutions
	int substitutions;

	//number of insertion
	int insertions;

	//number of deletion
	int deletions;

	//identity
	double identity;

} pytrf_ATR;

extern PyTypeObject pytrf_ATRType;

#endif
