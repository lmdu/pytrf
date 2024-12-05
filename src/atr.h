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

	//seed start position
	Py_ssize_t sstart;

	//seed end position
	Py_ssize_t send;

	//seed repeats
	int srepeat;

	//seed length
	int slen;

	//motif length
	int mlen;

	//tandem repeat
	float repeat;

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
	float identity;

} pytrf_ATR;

extern PyTypeObject pytrf_ATRType;

#endif
