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
#define MIN3(a, b, c) MIN(MIN((a), (b)), (c))
#define MIN4(a, b, c, d) MIN(MIN((a), (b)), MIN((c), (d)))

/*
 *@param s str, motif sequence
 *@param l int, motif length
 *@param m int, min motif length
 */
static int is_redundant_motif(char *s, int l, int m) {
	int i, j, b;

	if (m == 1) {
		return 0;
	}

	for (j = 1; j <= m; ++j) {
		b = l - j;
		i = 0;

		while ((i < b) && (s[i] == s[i+j])) {
			++i;
		}

		if (i == b) {
			return 1;
		}
	}

	return 0;
}

/*
 *@param s str, motif sequence
 *@param m int, motif length
 */
static void reverse_motif(char *s, int m) {
	int i, j;
	char b;

	for (i = 0, j = m - 1; i < j; ++i, --j) {
		b = s[i];
		s[i] = s[j];
		s[j] = b;
	}
}

/*
 * @param n int, maximum extend length
 * @param m int, maximum motif length
 * @return array, dynamic programming edit distance matrix
 */
static int** initial_matrix(int n, int m) {
	int i, j;
	int **d;

	d = (int **)malloc(sizeof(int *)*(n+1));

	for (i = 0; i <= n; ++i) {
		d[i] = (int *)malloc(sizeof(int)*(m+1));
		d[i][0] = i;
	}

	for (j = 0; j <= m; ++j) {
		d[0][j] = j;
	}

	return d;
}

void print_matrix(int **mx, int n, int m) {
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= m; ++j) {
			printf("%d\t", mx[i][j]);
		}
		printf("\n");
	}
}

/*
 * @param d array, edit distance matrix
 * @param n int, maximum extend length
 */
static void release_matrix(int **d, int n) {
	for (int i = 0; i <= n; ++i) {
		free(d[i]);
	}

	free(d);
}

/*
 * @param b char, current base in sequence
 * @param s str, motif sequence
 * @param m int, motif length
 * @param i int, row number
 * @param d array, wrap-around dynamic programming matrix
 * @return bool, 0 for no error occured, 1 for error occured
 */
static int wrap_around_distance(char b, char *s, int m, int i, int **d) {
	//column number
	int j;

	//match or mismatch cost, 0 or 1
	int c;

	//first pass
	//for the first column j = 1
	c = b == s[0] ? 0 : 1;
	d[i][1] = MIN3(d[i-1][0]+c, d[i-1][m]+c, d[i-1][1]+1);

	//for the column j > 1
	for (j = 2; j <= m; ++j) {
		c = b == s[j-1] ? 0 : 1;
		d[i][j] = MIN3(d[i-1][j-1]+c, d[i][j-1]+1, d[i-1][j]+1);
	}

	//sencond pass
	d[i][1] = MIN(d[i][1], d[i][m]+1);

	for (j = 2; j < m; ++j) {
		d[i][j] = MIN(d[i][j], d[i][j-1]+1);
	}

	//error occured or not
	return d[i][m] > d[i-1][m] ? 1 : 0;
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
 * @return int, extend length or row number of matrix
 */
static int wrap_around_extend(const char *s, char *ms, int ml, int **mx, Py_ssize_t st, int n, int me, int dr) {
	//row number of matrix
	int i;

	//consecutive errors
	int ce = 0;

	if (n <= 0) {
		return 0;
	}

	//fill the matrix row by row
	for (i = 1; i <= n; ++i) {
		if (wrap_around_distance(s[st+i*dr], ms, ml, i, mx) > 0) {
			++ce;
		} else {
			ce = 0;
		}

		if (ce > me) {
			break;
		}
	}

	if (i > n) {
		i = n;
	}

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
 * @param eds int array, number of mat, sub, del and ins
 * @return int, column number of minimum distance
 */
static float wrap_around_backtrace(int **mx, int m, int i, int dr, int *eds) {
	int mat = 0;
	int sub = 0;
	int del = 0;
	int ins = 0;

	//alignment ratio
	float r;

	//minimum value of edit distance
	int v;

	//column number of matrix
	int j;

	//extend length
	int l;

	j = m;
	l = i;

	if (i <= 0) {
		return 0;
	}

	while (i > 0 || j > 0) {
		//go back through second pass
		if ((i > 0) && (j > 0) && (j < m)) {
			if (j == 1) {
				if (mx[i][j] == (mx[i][m] + 1)) {
					
					++del;
					j = m;
					continue;
				}
			} else {
				if (mx[i][j] == (mx[i][j-1] + 1)) {
					++del;
					--j;
					continue;
				}
			}
		} else if (i == 0) {
			++del;
			--j;
			continue;
		}

		//go back through first pass
		if (j == 1) {
			v = MIN3(mx[i-1][m], mx[i-1][0], mx[i-1][1]);

			if (v == mx[i-1][m]) {
				if (v == mx[i][j]) {
					++mat;
				} else {
					++sub;
				}

				--i;
				j = m;
			} else if (v == mx[i-1][0]) {
				if (v == mx[i][j]) {
					++mat;
				} else {
					++sub;
				}

				--i;
				--j;
			} else if (v == (mx[i-1][1])) {
				++ins;
				--i;
			}
		} else {
			v = MIN3(mx[i-1][j-1], mx[i-1][j], mx[i][j-1]);

			if (v == mx[i-1][j-1]) {
				if (v == mx[i][j]) {
					++mat;
				} else {
					++sub;
				}

				--i;
				--j;
			} else if (v == mx[i-1][j]) {
				++ins;
				--i;
			} else if (v == mx[i][j-1]) {
				++del;
				--j;
			}
		}
	}

	eds[0] = mat;
	eds[1] = sub;
	eds[2] = ins;
	eds[3] = del;

	r = 1.0 * mat / l;

	return r;
}

static PyObject* pytrf_itrfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	int i;

	static char* keywords[] = {
		"chrom", "seq", "min_motif_size", "max_motif_size", "min_seed_repeat", "min_seed_length",
		"max_consecutive_error", "min_extend_identity", "max_extend_length", NULL
	};

	pytrf_ITRFinder *obj = (pytrf_ITRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	//initialize parameters
	obj->next_start = 0;
	obj->seed_minrep = 3;
	obj->seed_minlen = 10;
	obj->max_errors = 3;
	obj->min_motif = 1;
	obj->max_motif = 6;
	obj->min_identity = 0.7;
	obj->extend_maxlen = 1000;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiiiifi", keywords,
									&obj->seqname, &obj->seqobj, &obj->min_motif,
									&obj->max_motif, &obj->seed_minrep, &obj->seed_minlen,
									&obj->max_errors, &obj->min_identity, &obj->extend_maxlen))
	{
		return NULL;
	}

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);
	obj->motif = (char *)malloc(obj->max_motif + 1);
	obj->matrix = initial_matrix(obj->extend_maxlen, obj->max_motif);

	obj->boundary = (Py_ssize_t *)malloc(sizeof(Py_ssize_t) * (obj->max_motif+1));
	for (i = 0; i <= obj->max_motif; ++i) {
		obj->boundary[i] = obj->size - i;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

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

	if (self->boundary) {
		free(self->boundary);
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

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_length;
	int tandem_match;
	int tandem_substitute;
	int tandem_insert;
	int tandem_delete;
	float tandem_identity;
	
	int left_extend;
	int right_extend;

	float left_ratio;
	float right_ratio;

	int left_edits[4];
	int right_edits[4];

	Py_ssize_t b;

	for (i = self->next_start; i < self->size;  ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (j = self->min_motif; j <= self->max_motif; ++j) {
			b = self->boundary[j];

			while ((i < b) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			seed_length = i + j - seed_start;
			seed_repeat = seed_length / j;

			//corrected length
			seed_length = seed_repeat * j;

			if ((seed_repeat >= self->seed_minrep) && (seed_length >= self->seed_minlen)) {
				//get motif sequence
				memcpy(self->motif, self->seq + seed_start, j);
				self->motif[j] = '\0';

				if (is_redundant_motif(self->motif, j, self->min_motif)) {
					i = seed_start;
					continue;
				}

				seed_end = seed_start + seed_length - 1;

				//extend to left
				extend_start = seed_start;
				extend_maxlen = extend_start;

				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				//need to reverse motif before extend to left
				reverse_motif(self->motif, j);

				left_extend = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, -1);
				left_ratio = wrap_around_backtrace(self->matrix, j, left_extend, -1, left_edits);

				//need to recover motif after extend to left
				reverse_motif(self->motif, j);
				//tandem_start = extend_start - extend_len + 1;

				//extend to right
				extend_start = seed_end;
				extend_maxlen = self->size - extend_start - 1;
				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				right_extend = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, 1);
				right_ratio = wrap_around_backtrace(self->matrix, j, right_extend, 1, right_edits);


				//tandem_end = extend_start + extend_len + 1;

				//tandem_align = tandem_match + tandem_insert + tandem_substitute + tandem_delete;
				//tandem_identity = tandem_match * 1.0 / tandem_align;

				//tandem_match = seed_length;
				tandem_match = 0;
				tandem_substitute = 0;
				tandem_insert = 0;
				tandem_delete = 0;
				tandem_start = seed_start + 1;
				tandem_end = seed_end + 1;
				
				if (left_ratio >= self->min_identity) {
					tandem_match += left_edits[0];
					tandem_substitute += left_edits[1];
					tandem_insert += left_edits[2];
					tandem_delete += left_edits[3];
					tandem_start -= left_extend;
				}
				
				if (right_ratio >= self->min_identity) {
					tandem_match += right_edits[0];
					tandem_substitute += right_edits[1];
					tandem_insert += right_edits[2];
					tandem_delete += right_edits[3];
					tandem_end += right_extend;
				}

				tandem_length = tandem_end - tandem_start + 1;

				if (tandem_match) {
					tandem_identity = 1.0 * tandem_match / (tandem_length - seed_length) * 100;
				} else {
					tandem_identity = 0;
				}

				//create new atr element object
				pytrf_ATR *atr = PyObject_New(pytrf_ATR, &pytrf_ATRType);
				atr->motif = PyUnicode_FromString(self->motif);
				atr->mlen = j;
				atr->seqid = Py_NewRef(self->seqname);
				atr->seqobj = Py_NewRef(self->seqobj);
				atr->start = tandem_start;
				atr->end = tandem_end;
				atr->sstart = seed_start + 1;
				atr->send = seed_end + 1;
				atr->srepeat = seed_repeat;
				atr->repeat = (tandem_length - tandem_insert + tandem_delete) / j;
				atr->length = tandem_length;
				atr->matches = tandem_match;
				atr->substitutions = tandem_substitute;
				atr->insertions = tandem_insert;
				atr->deletions = tandem_delete;
				atr->identity = tandem_identity;

				self->next_start = tandem_end;

				return (PyObject *)atr;
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

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_length;
	int tandem_repeat;
	int tandem_match;
	int tandem_substitute;
	int tandem_insert;
	int tandem_delete;
	float tandem_identity;

	int left_extend;
	int right_extend;

	float left_ratio;
	float right_ratio;

	int left_edits[4];
	int right_edits[4];

	Py_ssize_t b;

	PyObject *itrs = PyList_New(0);
	PyObject *tmp;

	for (i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (j = self->min_motif; j <= self->max_motif; ++j) {
			b = self->boundary[j];

			while ((i < b) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			seed_length = i + j - seed_start;
			seed_repeat = seed_length / j;

			//corrected length
			seed_length = seed_repeat * j;

			if (seed_repeat >= self->seed_minrep && seed_length >= self->seed_minlen) {
				//get motif sequence
				memcpy(self->motif, self->seq + seed_start, j);
				self->motif[j] = '\0';

				if (is_redundant_motif(self->motif, j, self->min_motif)) {
					i = seed_start;
					continue;
				}

				seed_end = seed_start + seed_length - 1;

				//extend to left
				extend_start = seed_start;
				extend_maxlen = extend_start;

				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				reverse_motif(self->motif, j);

				left_extend = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, -1);

				left_ratio = wrap_around_backtrace(self->matrix, j, left_extend, -1, left_edits);

				reverse_motif(self->motif, j);

				//extend to right
				extend_start = seed_end;
				extend_maxlen = self->size - extend_start - 1;
				if (extend_maxlen > self->extend_maxlen) {
					extend_maxlen = self->extend_maxlen;
				}

				right_extend = wrap_around_extend(self->seq, self->motif, j, self->matrix, extend_start,
												extend_maxlen, self->max_errors, 1);

				right_ratio = wrap_around_backtrace(self->matrix, j, right_extend, 1, right_edits);

				tandem_match = 0;
				tandem_substitute = 0;
				tandem_insert = 0;
				tandem_delete = 0;
				tandem_start = seed_start + 1;
				tandem_end = seed_end + 1;

				if (left_ratio >= self->min_identity) {
					tandem_match += left_edits[0];
					tandem_substitute += left_edits[1];
					tandem_insert += left_edits[2];
					tandem_delete += left_edits[3];
					tandem_start -= left_extend;
				}
				
				if (right_ratio >= self->min_identity) {
					tandem_match += right_edits[0];
					tandem_substitute += right_edits[1];
					tandem_insert += right_edits[2];
					tandem_delete += right_edits[3];
					tandem_end += right_extend;
				}

				tandem_length = tandem_end - tandem_start + 1;
				tandem_repeat = (tandem_length - tandem_insert + tandem_delete) / j;

				if (tandem_match) {
					tandem_identity = 1.0 * tandem_match / (tandem_length - seed_length) * 100;
				} else {
					tandem_identity = 0;
				}

				tmp = Py_BuildValue("Onnsiinniiiiiif", self->seqname, seed_start + 1, seed_end + 1, self->motif,
									j, seed_repeat, tandem_start, tandem_end, tandem_repeat, tandem_length,
									tandem_match, tandem_substitute, tandem_insert, tandem_delete, tandem_identity);

				PyList_Append(itrs, tmp);
				Py_DECREF(tmp);

				i = tandem_end - 1;
				break;
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
