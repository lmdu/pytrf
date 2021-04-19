/*
 * itr.h -- imperfect tandem repeats
 *
 */

#define PY_SSIZE_T_CLEAN
#ifndef STRIPY_ITR_H
#define STRIPY_ITR_H
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

	//seed min repeats
	int seed_minrep;

	//seed min length
	int seed_minlen;

	//maximal continuous error
	int max_errors;

	//min motif size
	int min_motif;

	//max motif size
	int max_motif;

	//mismatch penalty
	//int mis_penalty;

	//indel penalty
	//int gap_penalty;

	//minimal score required to form a imperfect tandem repeat
	//int min_score;

	//maximal extend length
	int extend_maxlen;

} stripy_ITRMiner;

extern PyTypeObject stripy_ITRMinerType;

#endif
