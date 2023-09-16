/*
 * itr.c -- imperfect tandem repeat finder
 *
 */

#define PY_SSIZE_T_CLEAN
#include "atr.h"
#include "itr.h"
#include "math.h"
#include "compat.h"
#include "structmember.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MIN3(a, b, c) MIN(MIN(a, b), c)

void reverse_motif(char *ms, int ml) {
	int i, j;
	char b;

	for (i = 0, j = ml - 1; i < j; ++i, --j) {
		b = ms[i];
		ms[i] = ms[j];
		ms[j] = b;
	}
}

/*
 * @param n int, maximum extend length
 * @param m int, maximum motif length
 * @return array, dynamic programming matrix
 */
static int** initial_matrix(int n, int m) {
	int i;
	int j;
	int **mx;

	mx = (int **)malloc(sizeof(int *)*(n+1));

	for (i = 0; i <= n; ++i) {
		mx[i] = (int *)malloc(sizeof(int)*(m+1));
	}

	for (i = 0; i <= n; ++i) {
		mx[i][0] = i;
	}

	for (j = 0; j <= m; ++j) {
		mx[0][j] = j;
	}

	return mx;
}

/*void print_matrix(int **mx, int n, int m) {
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= m; ++j) {
			printf("%d\t", mx[i][j]);
		}
		printf("\n");
	}
}*/

/*
 * @param mx array, dp matrix
 * @param n int, maximum extend length
 */
static void release_matrix(int **mx, int n) {
	for (int i = 0; i <= n; ++i) {
		free(mx[i]);
	}

	free(mx);
}

/*
 * @param b char, current base in sequence
 * @param ms str, motif sequence
 * @param ml int, motif length
 * @param i int, row number
 * @param mx array, wrap-around dynamic programming matrix
 * @return int, column index of mimimum edit distance
 */
static int wrap_around_distance(char b, char *ms, int ml, int i, int **mx) {
	//column number
	int j;

	//match or mismatch cost, 0 or 1 
	int c;

	//column number of minimum edit distance
	int m;

	//first pass
	c = b == ms[0] ? 0 : 1;
	mx[i][1] = MIN3(mx[i-1][0]+c, mx[i-1][ml]+c, mx[i-1][1]+1);
	//mx[i][1] = MIN(mx[i-1][ml]+c, mx[i-1][1]+1);

	for (j = 2; j <= ml; ++j) {
		c = b == ms[j-1] ? 0 : 1;
		mx[i][j] = MIN3(mx[i-1][j-1]+c, mx[i][j-1]+1, mx[i-1][j]+1);
	}

	//sencond pass
	mx[i][1] = MIN(mx[i][1], mx[i][ml]+1);
	m = 1;

	for (j = 2; j < ml; ++j) {
		mx[i][j] = MIN(mx[i][j], mx[i][j-1]+1);

		if (mx[i][j] <= mx[i][j-1]) {
			m = j;
		}
	}

	return m;
}

/*
 * @param s str, DNA sequence
 * @param ms str, motif sequence
 * @param ml int, motif length
 * @param mx array, dynamic programming matrix
 * @param st int, start position to extend
 * @param n int, maximum allowed length to extend
 * @param me int, maximum allowed consecutive errors
 * @param dr int, extend direction -1 to left and 1 to right
 * @return int, extend length
 */
static int wrap_around_extend(const char *s, char *ms, int ml, int **mx, Py_ssize_t st, int n, int me, int dr) {
	//row number of matrix
	int i;

	//column number of minimum edit
	int j;

	//column number of minimum edit in upper row
	int k = 0;

	//consecutive errors
	int ce = 0;

	if (!n) {
		return 0;
	}

	for (i = 1; i <= n; ++i) {
		j = wrap_around_distance(s[st+i*dr], ms, ml, i, mx);
		ce = mx[i][j] > mx[i-1][k] ? ce + 1 : 0;

		if (ce > me) {
			break;
		}

		k = j;
	}

	i = i > n ? n : i;
	i -= ce;

	return i;
}

/*
 * @param s str, DNA sequence
 * @param ms str, motif sequence
 * @param ml int, motif length
 * @param mx array, dynamic programming matrix
 * @param st int, start position to extend
 * @param i int, row number in matrix
 * @param dr int, backtrace direction -1 to left and 1 to right
 * @param nm int, number of matches
 * @param ns int, number of substitutions
 * @param ni int, number of insertions
 * @param nd int, number of deletions
 * @return int, column number of minimum distance
 */
static int wrap_around_backtrace(const char *s, char *ms, int ml, int **mx, Py_ssize_t st, int i, int dr, int *nm, int *ns, int *ni, int *nd) {
	//match or mismatch cost
	int c;

	//column number of matrix
	int j;

	//column number for mimimum distance
	int k = 0;

	//find mimimum edit distance
	for (j = 1; j <= ml; ++j) {
		if (mx[i][j] <= mx[i][j-1]) {
			k = j;
		}
	}

	j = k;

	while (i > 0) {
		c = s[st+i*dr] == ms[j-1] ? 0 : 1;

		if (j == 1) {
			if (mx[i][1] == mx[i-1][ml]+c) {
				if (c == 0) {
					++*nm;
				} else {
					++*ns;
				}

				j = ml;
			} else if (mx[i][1] == mx[i][ml]+1) {
				++*nd;
				j = ml;
			} else if (mx[i][1] == mx[i-1][1]+1) {
				++*ni;
			}

			--i;
		} else {
			if (mx[i][j] == mx[i][j-1]+1) {
				++*nd;
				--j;
			} else if (mx[i][j] == mx[i-1][j-1]+c) {
				if (c == 0) {
					++*nm;
				} else {
					++*ns;
				}
				--i;
				--j;
			} else if (mx[i][j] == mx[i-1][j]+1) {
				++*ni;
				--i;
			}
		}
	}

	return k;
}

static PyObject* pytrf_itrfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {
		"chrom", "seq", "max_motif_size", "min_seed_repeat", "min_seed_length",
		"max_consecutive_error", "min_identity", "max_extend_length", NULL
	};

	pytrf_ITRFinder *obj = (pytrf_ITRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	//initialize parameters
	obj->next_start = 0;
	obj->seed_minrep = 3;
	obj->seed_minlen = 10;
	obj->max_errors = 3;
	obj->max_motif = 6;
	obj->min_identity = 70.0;
	obj->extend_maxlen = 1000;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiiidi", keywords,
									&obj->seqname, &obj->seqobj, &obj->max_motif,
									&obj->seed_minrep, &obj->seed_minlen, &obj->max_errors,
									&obj->min_identity, &obj->extend_maxlen))
	{
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);
	obj->motif = (char *)malloc(obj->max_motif + 1);
	obj->matrix = initial_matrix(obj->extend_maxlen, obj->max_motif);

	return (PyObject *)obj;
}

void pytrf_itrfinder_dealloc(pytrf_ITRFinder *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);

	self->seq = NULL;

	if (self->motif) {
		free(self->motif);
	}

	if (self->matrix) {
		release_matrix(self->matrix, self->extend_maxlen);
	}

	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* pytrf_itrfinder_repr(pytrf_ITRFinder *self) {
	return PyUnicode_FromFormat("<ATRFinder> for sequence %S", self->seqname);
}

static PyObject* pytrf_itrfinder_iter(pytrf_ITRFinder *self) {
	self->next_start = 0;

	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* pytrf_itrfinder_next(pytrf_ITRFinder *self) {
	int j;
	Py_ssize_t i;

	Py_ssize_t seed_start;
	Py_ssize_t seed_end;
	int seed_length;
	int seed_repeat;

	Py_ssize_t extend_start;
	int extend_maxlen;
	int extend_len;

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_align;
	int tandem_length;
	int tandem_match;
	int tandem_substitute = 0;
	int tandem_insert = 0;
	int tandem_delete = 0;
	double tandem_identity;

	Py_ssize_t boundary;

	for (i = self->next_start; i < self->size;  ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (j = 1; j <= self->max_motif; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			seed_length = i + j - seed_start;
			seed_repeat = seed_length/j;

			//corrected length
			seed_length = seed_repeat*j;

			if (seed_repeat >= self->seed_minrep && seed_length >= self->seed_minlen) {
				//get motif sequence
				memcpy(self->motif, self->seq + seed_start, j);
				self->motif[j] = '\0';

				seed_end = seed_start + seed_length - 1;
				tandem_match = seed_length;

				//extend to left
				extend_start = seed_start;
				extend_maxlen = extend_start;

				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				//need to reverse motif before extend to left
				reverse_motif(self->motif, j);

				extend_len = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, -1);

				if (extend_len > 0) {
					wrap_around_backtrace(self->seq, self->motif, j, self->matrix, extend_start,
											extend_len, -1, &tandem_match, &tandem_substitute,
											&tandem_insert, &tandem_delete);
				}
				//need to recover motif after extend to left
				reverse_motif(self->motif, j);
				tandem_start = extend_start - extend_len + 1;

				//extend to right
				extend_start = seed_end;
				extend_maxlen = self->size - extend_start - 1;
				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				extend_len = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, 1);

				if (extend_len > 0) {
					wrap_around_backtrace(self->seq, self->motif, j, self->matrix, extend_start,
											extend_len, 1, &tandem_match, &tandem_substitute,
											&tandem_insert, &tandem_delete);
				}

				tandem_end = extend_start + extend_len + 1;

				tandem_align = tandem_match + tandem_insert + tandem_substitute + tandem_delete;
				tandem_identity = (tandem_match * 1.0 / tandem_align) * 100;

				if (tandem_identity >= self->min_identity) {
					tandem_length = tandem_end - tandem_start + 1;

					//create new atr element object
					pytrf_ATR *atr = PyObject_New(pytrf_ATR, &pytrf_ATRType);
					atr->motif = PyUnicode_FromString(self->motif);
					atr->mlen = j;
					atr->seqid = Py_NewRef(self->seqname);
					atr->seqobj = Py_NewRef(self->seqobj);
					atr->start = tandem_start;
					atr->end = tandem_end;
					atr->length = tandem_length;
					atr->matches = tandem_match;
					atr->substitutions = tandem_substitute;
					atr->insertions = tandem_insert;
					atr->deletions = tandem_delete;
					atr->identity = tandem_identity;

					self->next_start = tandem_end;

					return (PyObject *)atr;
				}

				tandem_match = 0;
				tandem_substitute = 0;
				tandem_insert = 0;
				tandem_delete = 0;
			}
			i = seed_start;
		}
	}

	return NULL;
}
 
static PyObject* pytrf_itrfinder_as_list(pytrf_ITRFinder *self) {
	int j;
	Py_ssize_t i;
	
	Py_ssize_t seed_start;
	Py_ssize_t seed_end;
	int seed_length;
	int seed_repeat;

	Py_ssize_t extend_start;
	int extend_maxlen;
	int extend_len;

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_align;
	int tandem_length;
	int tandem_match;
	int tandem_substitute = 0;
	int tandem_insert = 0;
	int tandem_delete = 0;
	double tandem_identity;

	Py_ssize_t boundary;

	PyObject *itrs = PyList_New(0);
	PyObject *tmp;

	for (i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (j = 1; j <= self->max_motif; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			seed_length = i + j - seed_start;
			seed_repeat = seed_length/j;

			//corrected length
			seed_length = seed_repeat*j;

			if (seed_repeat >= self->seed_minrep && seed_length >= self->seed_minlen) {
				//get motif sequence
				memcpy(self->motif, self->seq + seed_start, j);
				self->motif[j] = '\0';

				seed_end = seed_start + seed_length - 1;
				tandem_match = seed_length;

				//extend to left
				extend_start = seed_start;
				extend_maxlen = extend_start;

				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				extend_len = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, -1);

				tandem_start = extend_start - extend_len + 1;
				if (extend_len > 0) {
					wrap_around_backtrace(self->seq, self->motif, j, self->matrix, extend_start,
											extend_len, -1, &tandem_match, &tandem_substitute,
											&tandem_insert, &tandem_delete);
				}

				//extend to right
				extend_start = seed_end;
				extend_maxlen = self->size - extend_start - 1;
				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				extend_len = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, 1);

				if (extend_len > 0) {
					wrap_around_backtrace(self->seq, self->motif, j, self->matrix, extend_start,
											extend_len, 1, &tandem_match, &tandem_substitute,
											&tandem_insert, &tandem_delete);
				}

				tandem_align = tandem_match + tandem_insert + tandem_substitute + tandem_delete;
				tandem_identity = (tandem_match * 1.0 / tandem_align) * 100;

				if (tandem_identity >= self->min_identity) {
					tandem_end = extend_start + extend_len + 1;
					tandem_length = tandem_end - tandem_start + 1;

					tmp = Py_BuildValue("Onnsiiiiiif", self->seqname, tandem_start, tandem_end, self->motif, j,
										tandem_length, tandem_match, tandem_substitute, tandem_insert,
										tandem_delete, tandem_identity);
					PyList_Append(itrs, tmp);
					Py_DECREF(tmp);

					i = tandem_end;
					break;
				}
				tandem_match = 0;
				tandem_substitute = 0;
				tandem_insert = 0;
				tandem_delete = 0;
			}
			i = seed_start;
		}
	}

	return itrs;
}

static PyMethodDef pytrf_itrfinder_methods[] = {
	{"as_list", (PyCFunction)pytrf_itrfinder_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject pytrf_ITRFinderType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "ATRFinder",
	.tp_basicsize = sizeof(pytrf_ITRFinder),
	.tp_dealloc = (destructor)pytrf_itrfinder_dealloc,
	.tp_repr = (reprfunc)pytrf_itrfinder_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "approximate tandem repeat finder",
	.tp_iter = (getiterfunc)pytrf_itrfinder_iter,
	.tp_iternext = (iternextfunc)pytrf_itrfinder_next,
	.tp_methods = pytrf_itrfinder_methods,
	.tp_new = pytrf_itrfinder_new,
};
