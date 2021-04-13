/* 
 * tre.h -- Tandem Repeat Element
 *
 * define the base class for tandem repeats
 * this class can be inherited by SSR and VNTR class
 *
 */

#ifndef STRIPY_TRE_H
#define STRIPY_TRE_H
#include "Python.h"

typedef struct {
	PyObject_HEAD

	//sequence name or identifier
	PyObject* seqid;

	//tandem repeat element motif sequence
	char* motif;

	//tandem repeat element motif length
	unsigned int mlen;

	//tandem repeat element start position
	Py_ssize_t start;

	//tandem repeat element stop position
	Py_ssize_t end;

	//number of tandem repeats
	unsigned int repeats;

	//tandem repeat element length
	unsigned int length;

} stripy_TRE;

extern PyTypeObject stripy_TREType;

#endif