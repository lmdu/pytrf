/*
 * itr.h -- imperfect tandem repeats
 *
 */

#define PY_SSIZE_T_CLEAN
#ifndef STRIA_ITR_H
#define STRIA_ITR_H
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

	//substitution penalty
	double sub_penalty;

	//insertion penalty
	double ins_penalty;

	//deletion penalty
	double del_penalty;

	//min match ratio
	double min_ratio;

	//maximal extend length
	int extend_maxlen;

} stria_ITRMiner;

typedef struct {
	PyObject_HEAD

	//input sequence name
	PyObject *seqid;

	//imperfect tandem repeat start position
	Py_ssize_t start;

	//imperfect tandem repeat stop position
	Py_ssize_t end;

	//motif sequence
	char *motif;

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

} stria_ITR;

extern PyTypeObject stria_ITRMinerType;
extern PyTypeObject stria_ITRType;

#endif
