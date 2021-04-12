/* 
 * tre.h -- Tandem Repeat Element
 *
 * define the base class for tandem repeats
 * this class can be inherited by SSR and VNTR class
 *
 */

#ifndef STRIPY_TRE_H
#define STRIPY_TRE_H

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

} stripy_TRE;

extern PyTypeObject stripy_TREType;

#endif