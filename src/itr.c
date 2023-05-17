/*
 * itr.c -- imperfect tandem repeat finder
 *
 */

#define PY_SSIZE_T_CLEAN
#include "atr.h"
#include "itr.h"
#include "math.h"
#include "structmember.h"

static int min_three(int a, int b, int c) {
	int d;
	d = a < b ? a : b;
	return d < c ? d : c;
}

static int min_two(int a, int b) {
	return a > b ? a : b;
}

static int** initial_matrix(int size, int mlen) {
	int i;
	int j;
	int **matrix;

	matrix = (int **)malloc(sizeof(int *)*size);

	for (i = 0; i <= size; ++i) {
		matrix[i] = (int *)malloc(sizeof(int)*mlen);
	}

	for (j = 0; j <= mlen; ++j) {
		matrix[0][j] = j;
	}

	for (i = 0; i <= size; ++i) {
		matrix[i][0] = i;
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

static int wrap_around_edit(char base, char *motif, int mlen, int row, int **matrix) {
	int col;
	int cost;
	int min_edits;

	//first pass
	cost = base == motif[0] ? 0 : 1;
	matrix[row][1] = min_three(matrix[row-1][0]+cost, matrix[row-1][mlen]+cost, matrix[row-1][1]+1);

	for (col = 2; col <= mlen; ++col) {
		cost = base == motif[col-1] ? 0 : 1;
		matrix[row][col] = min_three(matrix[row-1][col-1]+cost, matrix[row][col-1]+1, matrix[row-1][col]+1);
	}

	//sencond pass
	matrix[row][1] = min_two(matrix[row][1], matrix[row][mlen]+1);
	min_edits = matrix[row][1];

	for (col = 2; col < mlen; ++col) {
		matrix[row][col] = min_two(matrix[row][col], matrix[row][col-1]+1);

		if (matrix[row][col] < cur_edits) {
			min_edits = matrix[row][col];
		}
	}

	return min_edits;
}

static int extend_to_left(const char *seq, char *motif, int mlen, int **matrix, Py_ssize_t start, int size, int max_error) {
	int i, j;

	//upper and current edits
	int upper_edits = 0;
	int cur_edits = 0;

	//successive errors
	int error = 0;

	char base;

	for (i = 1; i <= size; ++i) {
		base = seq[start-i];
		cur_edits = wrap_around_edit(base, motif, mlen, i, matrix);

		if (cur_edits > upper_edits) {
			++error;
		} else {
			error = 0;
		}

		if (error > max_error) {
			break;
		}

		upper_edits = cur_edits;
	}

	return i - error;
}

static int extend_to_right(const char *seq, char *motif, int mlen, int **matrix, Py_ssize_t start, int size, int max_error) {
	int i, j;

	//upper and current edits
	int upper_edits = 0;
	int cur_edits = 0;

	//successive errors
	int error = 0;

	char base;

	for (i = 1; i <= size; ++i) {
		base = seq[start+i];
		cur_edits = wrap_around_edit(base, motif, mlen, i, matrix);

		if (cur_edits > upper_edits) {
			++error;
		} else {
			error = 0;
		}

		if (error > max_error) {
			break;
		}

		upper_edits = cur_edits;
	}

	return i - error;
}

static int wrap_around_backtrace(int **matrix, int row, int mlen, int *mat, int *sub, int *ins, int *del) {
	
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

static PyObject* pytrf_itrfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {
		"name", "seq", "min_motif_size", "max_motif_size", "seed_min_repeat",
		"seed_min_length", "max_continuous_errors", "substitution_penalty",
		"insertion_penalty", "deletion_penalty", "min_match_ratio", "max_extend_length",
		NULL
	};

	pytrf_ITRFinder *obj = (pytrf_ITRFinder *)type->tp_alloc(type, 0);
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
	obj->motif = NULL;
	obj->matrix = NULL;

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
	return PyUnicode_FromFormat("<ITRFinder> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* pytrf_itrfinder_iter(pytrf_ITRFinder *self) {
	self->next_start = 0;

	if (!self->motif) {
		self->motif = (char *)malloc(self->max_motif + 1);
	}

	if (!self->matrix) {
		self->matrix = initial_matrix(self->extend_maxlen);
	}

	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* pytrf_itrfinder_next(pytrf_ITRFinder *self) {
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

	Py_ssize_t boundary;

	for (Py_ssize_t i = self->next_start; i < self->size;  ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
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
					memcpy(self->motif, self->seq + seed_start, j);
					self->motif[j] = '\0';

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

					extend_end = build_left_matrix(self->seq, self->motif, j, self->matrix, extend_start,
												   extend_maxlen, self->max_errors);
					extend_len = backtrace_matrix(self->matrix, extend_end, &tandem_match, &substitution,
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

						extend_end = build_right_matrix(self->seq, self->motif, j, self->matrix, extend_start,
														extend_maxlen, self->max_errors);
						extend_len = backtrace_matrix(self->matrix, extend_end, &tandem_match, &substitution,
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

							//create new atr element object
							pytrf_ATR *atr = PyObject_New(pytrf_ATR, &pytrf_ATRType);
							atr->motif = (char *)malloc(j + 1);
							memcpy(atr->motif, self->motif, j);
							atr->motif[j] = '\0';
							atr->mlen = j;
							atr->seqid = self->seqname;
							Py_INCREF(atr->seqid);
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
					}
				}
			}
			i = seed_start;
		}
	}

	free(self->motif);
	release_matrix(self->matrix, self->extend_maxlen);
	self->motif = NULL;
	self->matrix = NULL;
	return NULL;
}
 
static PyObject* pytrf_itrfinder_as_list(pytrf_ITRFinder *self) {
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

	if (!self->motif) {
		self->motif = (char *)malloc(self->max_motif + 1);
	}

	if (!self->matrix) {
		self->matrix = initial_matrix(self->extend_maxlen);
	}

	Py_ssize_t boundary;

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		seed_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
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
					memcpy(self->motif, self->seq + seed_start, j);
					self->motif[j] = '\0';

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

					extend_end = build_left_matrix(self->seq, self->motif, j, self->matrix, extend_start,
												   extend_maxlen, self->max_errors);
					extend_len = backtrace_matrix(self->matrix, extend_end, &tandem_match, &substitution,
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

						extend_end = build_right_matrix(self->seq, self->motif, j, self->matrix, extend_start,
														extend_maxlen, self->max_errors);
						extend_len = backtrace_matrix(self->matrix, extend_end, &tandem_match, &substitution,
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

							tmp = Py_BuildValue("Onnsiiiiiif", self->seqname, tandem_start, tandem_end, self->motif, j,
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

	free(self->motif);
	release_matrix(self->matrix, self->extend_maxlen);
	return itrs;
}

static PyMethodDef pytrf_itrfinder_methods[] = {
	{"as_list", (PyCFunction)pytrf_itrfinder_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject pytrf_ITRFinderType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "ITRFinder",
	.tp_basicsize = sizeof(pytrf_ITRFinder),
	.tp_dealloc = (destructor)pytrf_itrfinder_dealloc,
	.tp_repr = (reprfunc)pytrf_itrfinder_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "imperfect or approximate tandem repeat finder",
	.tp_iter = (getiterfunc)pytrf_itrfinder_iter,
	.tp_iternext = (iternextfunc)pytrf_itrfinder_next,
	.tp_methods = pytrf_itrfinder_methods,
	.tp_new = pytrf_itrfinder_new,
};
