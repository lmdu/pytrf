/*
 * itr.c -- imperfect tandem repeat
 *
 */

#define PY_SSIZE_T_CLEAN
#include "itr.h"
#include "math.h"
#include "structmember.h"

static int min_three(int a, int b, int c) {
	int d;
	d = a < b ? a : b;
	return d < c ? d : c;
}

static int** initial_matrix(int size) {
	int i;
	int **matrix = (int **)malloc(sizeof(int *)*size);

	for (i = 0; i <= size; ++i) {
		matrix[i] = (int *)malloc(sizeof(int)*size);
	}

	for (i = 0; i <= size; ++i) {
		matrix[i][0] = i;
		matrix[0][i] = i;
	}

	return matrix;
}

static void release_matrix(int **matrix, int size) {
	for (int i = 0; i <= size; ++i) {
		free(matrix[i]);
	}
	free(matrix);
}

static int* build_left_matrix(const char *seq, char *motif, int mlen, int **matrix, Py_ssize_t start, int size, int max_error) {
	char h; //horizantal base in matrix
	char v; //vertical base in matrix
	int i = 0;
	int j = 0;
	int x = 0;
	int y = 0;
	int last_x = 0;
	int last_y = 0;
	int error = 0; //consective errors
	int smaller;
	
	static int res[2]; //result arrary

	for (x = 1, y = 1; x <= size && y <= size; ++x, ++y) {
		h = seq[start-y];
		v = motif[(mlen-x%mlen)%mlen];
		
		//fill column, column number fixed
		if (i != y) {
			for (i = 1; i < x; ++i){
				if (h == motif[(mlen-i%mlen)%mlen]) {
					matrix[i][y] = matrix[i-1][y-1];
				} else {
					matrix[i][y] = min_three(matrix[i-1][y-1], matrix[i-1][y], matrix[i][y-1]) + 1;
				}
			}
		}
		//fill row, row number fixed
		if (j != x) {
			for (j = 1; j < y; ++j) {
				if (v == seq[start-j]) {
					matrix[x][j] = matrix[x-1][j-1];
				} else {
					matrix[x][j] = min_three(matrix[x-1][j-1], matrix[x-1][j], matrix[x][j-1]) + 1;
				}
			}
		}

		i = y;
		j = x;

		if (h == v) {
			matrix[x][y] = matrix[x-1][y-1];
			error = 0;
		} else {
			if (error == 0) {
				last_x = x - 1;
				last_y = y - 1;
			}
			
			error++;

			if (error > max_error) {
				break;
			}
			
			matrix[x][y] = min_three(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) + 1;
		}
		
		smaller = min_three(matrix[x][y], matrix[x-1][y], matrix[x][y-1]);
		if (smaller != matrix[x][y]) {
			if (matrix[x-1][y] != matrix[x][y-1]) {
				if (smaller == matrix[x][y-1]) {
					y -= 1;
				} else {
					x -= 1;
				}
			}
		}
	}

	if (error) {
		res[0] = last_x;
		res[1] = last_y;
	} else {
		res[0] = --x;
		res[1] = --y;
	}
	return res;
}

static int* build_right_matrix(const char *seq, char *motif, int mlen, int **matrix, Py_ssize_t start, int size, int max_error) {
	char h; //horizantal base in sequence
	char v; //vertical base
	int i = 0;
	int j = 0;
	int x = 0;
	int y = 0;
	int last_x = 0;
	int last_y = 0;
	int error = 0; //consective errors
	int smaller;
	
	static int res[2]; //result arrary

	for (x=1, y=1; x <= size && y <= size; ++x, ++y) {
		h = seq[start+y];
		v = motif[(x-1)%mlen];
		
		//fill column, column number fixed
		if (i != y) {
			for (i = 1; i < x; ++i) {
				if (h == motif[(i-1)%mlen]) {
					matrix[i][y] = matrix[i-1][y-1];
				} else {
					matrix[i][y] = min_three(matrix[i-1][y-1], matrix[i-1][y], matrix[i][y-1]) + 1;
				}
			}
		}
		//fill row, row number fixed
		if (j != x) {
			for (j = 1; j < y; ++j) {
				if (v == seq[start+j]) {
					matrix[x][j] = matrix[x-1][j-1];
				} else {
					matrix[x][j] = min_three(matrix[x-1][j-1], matrix[x-1][j], matrix[x][j-1]) + 1;
				}
			}
		}

		i = y;
		j = x;

		if (h == v) {
			matrix[x][y] = matrix[x-1][y-1];
			error = 0;
		} else {
			if (error == 0) {
				last_x = x - 1;
				last_y = y - 1;
			}
			
			error++;

			if (error > max_error) {
				break;
			}
			
			matrix[x][y] = min_three(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) + 1;
		}

		smaller = min_three(matrix[x][y], matrix[x-1][y], matrix[x][y-1]);
		if (smaller != matrix[x][y]) {
			if (matrix[x-1][y] != matrix[x][y-1]) {
				if (smaller == matrix[x][y-1]) {
					y -= 1;
				} else {
					x -= 1;
				}
			}
		}
	}

	if (error) {
		res[0] = last_x;
		res[1] = last_y;
	} else {
		res[0] = --x;
		res[1] = --y;
	}
	return res;
}

static int backtrace_matrix(int **matrix, int *diagonal, int *mat, int *sub, int *ins, int *del) {
	int i = *diagonal;
	int j = *(diagonal+1);
	int cost;
	int r = j;

	while(i > 0 && j > 0) {
		cost = min_three(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]);
		if (cost == matrix[i-1][j-1]) {
			if (cost == matrix[i][j]) {
				*mat += 1;
			} else {
				*sub += 1;
			}
			i--;
			j--;
		} else if (cost == matrix[i-1][j]){
			*del += 1;
			i--;
		} else {
			*ins += 1;
			j--;
		}
	}

	if (i>0) {
		*del += 1;
	} else if (j>0) {
		*ins += 1;
	}

	return r;
}

static PyObject* stria_itrminer_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {
		"name", "seq", "min_motif_size", "max_motif_size", "seed_min_repeat",
		"seed_min_length", "max_continuous_errors", "substitution_penalty",
		"insertion_penalty", "deletion_penalty", "min_match_ratio", "max_extend_length",
		NULL
	};

	stria_ITRMiner *obj = (stria_ITRMiner *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	//initialize parameters
	obj->next_start = 0;
	obj->seed_minrep = 3;
	obj->seed_minlen = 10;
	obj->max_errors = 2;
	obj->min_motif = 1;
	obj->max_motif = 6;
	obj->sub_penalty = 0.5;
	obj->ins_penalty = 1.0;
	obj->del_penalty = 1.0;
	obj->min_ratio = 0.7;
	obj->extend_maxlen = 2000;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiiiiddddi", keywords,
									&obj->seqname, &obj->seqobj, &obj->min_motif,
									&obj->max_motif, &obj->seed_minrep, &obj->seed_minlen,
									&obj->max_errors, &obj->sub_penalty, &obj->ins_penalty,
									&obj->del_penalty, &obj->min_ratio, &obj->extend_maxlen))
	{
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	return (PyObject *)obj;
}

void stria_itrminer_dealloc(stria_ITRMiner *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* stria_itrminer_repr(stria_ITRMiner *self) {
	return PyUnicode_FromFormat("<ITRMiner> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* stria_itrminer_iter(stria_ITRMiner *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* stria_itrminer_next(stria_ITRMiner *self) {
	Py_ssize_t seed_start;
	Py_ssize_t seed_end;
	int seed_length;
	int seed_repeat;
	int seed_good;

	//int matches;
	int substitution;
	int insertion;
	int deletion;

	Py_ssize_t extend_start;
	int* extend_end;
	int extend_maxlen;
	int extend_len;

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_length;
	int tandem_match;
	int tandem_substitute;
	int tandem_insert;
	int tandem_delete;
	double tandem_identity;

	double align_rate;

	char* motif = (char *)malloc(self->max_motif + 1);

	int **matrix = initial_matrix(self->extend_maxlen);

	for (Py_ssize_t i = self->next_start; i <= self->size;  ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			seed_length = i + j - seed_start;
			seed_repeat = seed_length/j;

			//corrected length
			seed_length = seed_repeat*j;

			if (seed_repeat >= self->seed_minrep && seed_length >= self->seed_minlen) {
				const char *p = self->seq + seed_start;
				seed_good = 1;

				for (int k = 1; k < self->min_motif; ++k) {
					int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l+k == j) {
						seed_good = 0;
						break;
					}
				}

				if (seed_good) {
					//get motif sequence
					memcpy(motif, self->seq + seed_start, j);
					motif[j] = '\0';

					seed_end = seed_start + seed_length - 1;
					tandem_match = seed_length;
					insertion = 0;
					deletion = 0;
					substitution = 0;

					//extend left flank
					extend_start = seed_start;
					extend_maxlen = extend_start;

					if (extend_maxlen > self->extend_maxlen) {
						extend_maxlen = self->extend_maxlen;
					}

					extend_end = build_left_matrix(self->seq, motif, j, matrix, extend_start,
												   extend_maxlen, self->max_errors);
					extend_len = backtrace_matrix(matrix, extend_end, &tandem_match, &substitution,
												  &insertion, &deletion);
					
					if (extend_len > 0) {
						align_rate = 1 - (substitution*self->sub_penalty + insertion*self->ins_penalty
									 + deletion*self->del_penalty) / extend_len;
					} else {
						align_rate = 1;
					}

					//if left is ok, extend to right
					if (align_rate >= self->min_ratio) {
						tandem_start = extend_start - extend_len + 1;
						tandem_substitute = substitution;
						tandem_insert = insertion;
						tandem_delete = deletion;

						substitution = 0;
						insertion = 0;
						deletion = 0;

						//extend right flank
						extend_start = seed_end;
						extend_maxlen = self->size - extend_start - 1;
						if (extend_maxlen > self->extend_maxlen) {
							extend_maxlen = self->extend_maxlen;
						}

						extend_end = build_right_matrix(self->seq, motif, j, matrix, extend_start,
														extend_maxlen, self->max_errors);
						extend_len = backtrace_matrix(matrix, extend_end, &tandem_match, &substitution,
													  &insertion, &deletion);
						
						//calcuate the alignment rate of extended sequence
						if (extend_len > 0) {
							align_rate = 1 - (substitution*self->sub_penalty + insertion*self->ins_penalty
										 + deletion*self->del_penalty) / extend_len;
						} else {
							align_rate = 1;
						}

						if (align_rate >= self->min_ratio) {
							tandem_end = extend_start + extend_len + 1;
							tandem_length = tandem_end - tandem_start + 1;
							tandem_substitute += substitution;
							tandem_insert += insertion;
							tandem_delete += deletion;
							tandem_identity = (tandem_match * 1.0 / tandem_length)*100;

							//create new itr element object
							stria_ITR *itr = PyObject_New(stria_ITR, &stria_ITRType);
							itr->motif = (char *)malloc(j + 1);
							memcpy(itr->motif, motif, j);
							itr->motif[j] = '\0';
							itr->mlen = j;
							itr->seqid = self->seqname;
							Py_INCREF(itr->seqid);
							itr->start = tandem_start;
							itr->end = tandem_end;
							itr->length = tandem_length;
							itr->matches = tandem_match;
							itr->substitutions = tandem_substitute;
							itr->insertions = tandem_insert;
							itr->deletions = tandem_delete;
							itr->identity = tandem_identity;

							self->next_start = tandem_end;
							return (PyObject *)itr;
						}
					}
				}
			}
			i = seed_start;
		}
	}

	free(motif);
	release_matrix(matrix, self->extend_maxlen);
	return NULL;
}
 
static PyObject* stria_itrminer_as_list(stria_ITRMiner *self) {
	PyObject *itrs = PyList_New(0);
	PyObject *tmp;
	
	Py_ssize_t seed_start;
	Py_ssize_t seed_end;
	int seed_length;
	int seed_repeat;
	int seed_good;

	//int matches;
	int substitution;
	int insertion;
	int deletion;

	Py_ssize_t extend_start;
	int* extend_end;
	int extend_maxlen;
	int extend_len;

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_length;
	int tandem_match;
	int tandem_substitute;
	int tandem_insert;
	int tandem_delete;
	double tandem_identity;

	double align_rate;
	
	char* motif = (char *)malloc(self->max_motif + 1);

	int **matrix = initial_matrix(self->extend_maxlen);

	for (Py_ssize_t i = 0; i <= self->size;  ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			seed_length = i + j - seed_start;
			seed_repeat = seed_length/j;

			//corrected length
			seed_length = seed_repeat*j;

			if (seed_repeat >= self->seed_minrep && seed_length >= self->seed_minlen) {
				const char *p = self->seq + seed_start;
				seed_good = 1;

				for (int k = 1; k < self->min_motif; ++k) {
					int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l+k == j) {
						seed_good = 0;
						break;
					}
				}

				if (seed_good) {
					//get motif sequence
					memcpy(motif, self->seq + seed_start, j);
					motif[j] = '\0';

					seed_end = seed_start + seed_length - 1;
					tandem_match = seed_length;
					insertion = 0;
					deletion = 0;
					substitution = 0;

					//extend left flank
					extend_start = seed_start;
					extend_maxlen = extend_start;

					if (extend_maxlen > self->extend_maxlen) {
						extend_maxlen = self->extend_maxlen;
					}

					extend_end = build_left_matrix(self->seq, motif, j, matrix, extend_start,
												   extend_maxlen, self->max_errors);
					extend_len = backtrace_matrix(matrix, extend_end, &tandem_match, &substitution,
												  &insertion, &deletion);
					
					if (extend_len > 0) {
						align_rate = 1 - (substitution*self->sub_penalty + insertion*self->ins_penalty
									 + deletion*self->del_penalty) / extend_len;
					} else {
						align_rate = 1;
					}

					//if left is ok, extend to right
					if (align_rate >= self->min_ratio) {
						tandem_start = extend_start - extend_len + 1;
						tandem_substitute = substitution;
						tandem_insert = insertion;
						tandem_delete = deletion;

						substitution = 0;
						insertion = 0;
						deletion = 0;

						//extend right flank
						extend_start = seed_end;
						extend_maxlen = self->size - extend_start - 1;
						if (extend_maxlen > self->extend_maxlen) {
							extend_maxlen = self->extend_maxlen;
						}

						extend_end = build_right_matrix(self->seq, motif, j, matrix, extend_start,
														extend_maxlen, self->max_errors);
						extend_len = backtrace_matrix(matrix, extend_end, &tandem_match, &substitution,
													  &insertion, &deletion);
						
						//calcuate the alignment rate of extended sequence
						if (extend_len > 0) {
							align_rate = 1 - (substitution*self->sub_penalty + insertion*self->ins_penalty
										 + deletion*self->del_penalty) / extend_len;
						} else {
							align_rate = 1;
						}

						if (align_rate >= self->min_ratio) {
							tandem_end = extend_start + extend_len + 1;
							tandem_length = tandem_end - tandem_start + 1;
							tandem_substitute += substitution;
							tandem_insert += insertion;
							tandem_delete += deletion;
							tandem_identity = (tandem_match * 1.0 / tandem_length)*100;

							tmp = Py_BuildValue("Onnsiiiiiif", self->seqname, tandem_start, tandem_end, motif, j,
												tandem_length, tandem_match, tandem_substitute, tandem_insert,
												tandem_delete, tandem_identity);
							PyList_Append(itrs, tmp);
							Py_DECREF(tmp);

							i = tandem_end;
							break;
						}
					}
				}
			}
			i = seed_start;
		}
	}

	free(motif);
	release_matrix(matrix, self->extend_maxlen);
	return itrs;
}

static PyMethodDef stria_itrminer_methods[] = {
	{"as_list", (PyCFunction)stria_itrminer_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject stria_ITRMinerType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"ITRMiner",                        /* tp_name */
	sizeof(stria_ITRMiner),          /* tp_basicsize */
	0,                              /* tp_itemsize */
	(destructor)stria_itrminer_dealloc,   /* tp_dealloc */
	0,                              /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_reserved */
	(reprfunc)stria_itrminer_repr,                              /* tp_repr */
	0,                              /* tp_as_number */
	0,                   /* tp_as_sequence */
	0,                   /* tp_as_mapping */
	0,                              /* tp_hash */
	0,                              /* tp_call */
	0,                              /* tp_str */
	0,                              /* tp_getattro */
	0,                              /* tp_setattro */
	0,                              /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,             /* tp_flags */
	"find microsatellites from DNA sequence",                              /* tp_doc */
	0,                              /* tp_traverse */
	0,                              /* tp_clear */
	0,                              /* tp_richcompare */
	0,                              /* tp_weaklistoffset */
	(getiterfunc)stria_itrminer_iter,     /* tp_iter */
	(iternextfunc)stria_itrminer_next,    /* tp_iternext */
	stria_itrminer_methods,          /* tp_methods */
	0,          /* tp_members */
	0,                               /* tp_getset */
	0,                              /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	0,                              /* tp_init */
	0,            /* tp_alloc */
	stria_itrminer_new,              /* tp_new */
};

/*
 * ITR element
 */

void stria_itr_dealloc(stria_ITR *self) {
	free(self->motif);
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* stria_itr_repr(stria_ITR *self) {
	return PyUnicode_FromFormat("<ITR> %s @ %s:%zd-%zd", self->motif,
								PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}

PyObject* stria_itr_as_list(stria_ITR *self) {
	return Py_BuildValue("Onnsiiiiiif", self->seqid, self->start, self->end, self->motif, self->mlen,
						self->length, self->matches, self->substitutions, self->insertions,
						self->deletions, self->identity);
}

PyObject* stria_itr_as_dict(stria_ITR *self) {
	return Py_BuildValue("{s:O,s:n,s:n,s:s,s:i,s:i,s:i,s:i,s:i,s:i,s:d}", "chrom", self->seqid,
						 "start", self->start, "end", self->end, "motif", self->motif, "type", self->mlen,
						 "length", self->length, "matches", self->matches, "substitutions", self->substitutions,
						 "insertions", self->insertions, "deletions", self->deletions, "identity", self->identity);
}

PyObject* stria_itr_as_gff(stria_ITR *self) {
	return PyUnicode_FromFormat("%S\tstria\tITR\t%zd\t%zd\t.\t+\t.\tMotif=%s;Type=%d;Length=%d;Match=%d;"
								"Substitutions=%d;Insertions=%d;Deletions=%dIdentity=%R\n", self->seqid, self->start,
								self->end, self->motif, self->mlen, self->length, self->matches, self->substitutions, 
								self->insertions, self->deletions, PyFloat_FromDouble(self->identity));
}

PyObject* stria_itr_as_string(stria_ITR *self, PyObject *args, PyObject *kwargs) {
	char *separator = "\t";
	char *terminator = "";
	static char* keywords[] = {"separator", "terminator", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|ss", keywords, &separator, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%S%s%zd%s%zd%s%s%s%d%s%d%s%d%s%d%s%d%s%d%s%R%s", self->seqid, separator,
								self->start, separator, self->end, separator, self->motif, separator,
								self->mlen, separator, self->length, separator, self->matches, separator,
								self->substitutions, separator, self->insertions, separator, self->deletions,
								separator, PyFloat_FromDouble(self->identity), terminator);
}

PyObject *stria_itr_get_seq(stria_ITR *self, void* closure) {
	PyObject* ret = PyUnicode_New(self->length, 127);
	Py_UCS1* p = PyUnicode_1BYTE_DATA(ret);
	memcpy(p, PyUnicode_AsUTF8(self->seqid)+self->start-1, self->length);
	return ret;
}

static PyMethodDef stria_itr_methods[] = {
	{"as_list", (PyCFunction)stria_itr_as_list, METH_NOARGS, NULL},
	{"as_dict", (PyCFunction)stria_itr_as_dict, METH_NOARGS, NULL},
	{"as_gff", (PyCFunction)stria_itr_as_gff, METH_NOARGS, NULL},
	{"as_string", (PyCFunction)stria_itr_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef stria_itr_getsets[] = {
	{"seq", (getter)stria_itr_get_seq, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef stria_itr_members[] = {
	{"chrom", T_OBJECT, offsetof(stria_ITR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(stria_ITR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(stria_ITR, end), READONLY},
	{"motif", T_STRING, offsetof(stria_ITR, motif), READONLY},
	{"type", T_INT, offsetof(stria_ITR, mlen), READONLY},
	{"length", T_INT, offsetof(stria_ITR, length), READONLY},
	{"matches", T_INT, offsetof(stria_ITR, matches), READONLY},
	{"substitutions", T_INT, offsetof(stria_ITR, substitutions), READONLY},
	{"insertions", T_INT, offsetof(stria_ITR, insertions), READONLY},
	{"deletions", T_INT, offsetof(stria_ITR, deletions), READONLY},
	{"identity", T_DOUBLE, offsetof(stria_ITR, identity), READONLY},
	{NULL}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
};

PyTypeObject stria_ITRType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"IRE",                        /* tp_name */
	sizeof(stria_ITR),          /* tp_basicsize */
	0,                              /* tp_itemsize */
	(destructor)stria_itr_dealloc,   /* tp_dealloc */
	0,                              /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_reserved */
	(reprfunc)stria_itr_repr,                              /* tp_repr */
	0,                              /* tp_as_number */
	0,                   /* tp_as_sequence */
	0,                   /* tp_as_mapping */
	0,                              /* tp_hash */
	0,                              /* tp_call */
	0,                              /* tp_str */
	0,                              /* tp_getattro */
	0,                              /* tp_setattro */
	0,                              /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,             /* tp_flags */
	"tandem repeat element",                              /* tp_doc */
	0,                              /* tp_traverse */
	0,                              /* tp_clear */
	0,                              /* tp_richcompare */
	0,                              /* tp_weaklistoffset */
	0,     /* tp_iter */
	0,    /* tp_iternext */
	stria_itr_methods,          /* tp_methods */
	stria_itr_members,          /* tp_members */
	stria_itr_getsets,                               /* tp_getset */
	0,                              /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	0,                              /* tp_init */
	0,            /* tp_alloc */
	PyType_GenericNew,              /* tp_new */
};
