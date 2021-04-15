/*
 * itr.c -- imperfect tandem repeat
 *
 */

#define PY_SSIZE_T_CLEAN
#include "itr.h"
#include "structmember.h"

static int min_three(int a, int b, int c){
	int d;
	d = a<b?a:b;
	return d<c?d:c;
}

static int** initial_matrix(unsigned int size){
	int i;
	int **matrix = (int **)malloc(sizeof(int *)*size);

	for(i=0; i<=size; i++){
		matrix[i] = (int *)malloc(sizeof(int)*size);
	}

	for(i=0; i<=size; ++i){
		matrix[i][0] = i;
		matrix[0][i] = i;
	}

	return matrix;
}

static void release_matrix(int **matrix, unsigned int size){
	for(int i=0; i<=size; ++i){
		free(matrix[i]);
	}
	free(matrix);
}

static int* build_left_matrix(const char *seq, char *motif, int mlen, int **matrix, Py_ssize_t start, int size, int max_error){
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

	for(x = 1, y = 1; x <= size && y <= size; ++x, ++y){
		h = seq[start-y];
		v = motif[(mlen-x%mlen)%mlen];
		
		//fill column, column number fixed
		if(i != y){
			for(i=1; i<x; i++){
				if(h == motif[(mlen-i%mlen)%mlen]){
					matrix[i][y] = matrix[i-1][y-1];
				}else{
					matrix[i][y] = min_three(matrix[i-1][y-1], matrix[i-1][y], matrix[i][y-1]) + 1;
				}
			}
		}
		//fill row, row number fixed
		if(j != x){
			for(j=1; j<y; j++){
				if(v == seq[start-j]){
					matrix[x][j] = matrix[x-1][j-1];
				}else{
					matrix[x][j] = min_three(matrix[x-1][j-1], matrix[x-1][j], matrix[x][j-1]) + 1;
				}
			}
		}

		i = y;
		j = x;

		if(h == v){
			matrix[x][y] = matrix[x-1][y-1];
			error = 0;
		}else{
			if(error == 0){
				last_x = x - 1;
				last_y = y - 1;
			}
			
			error++;

			if(error > max_error){
				break;
			}
			
			matrix[x][y] = min_three(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) + 1;
		}
		
		smaller = min_three(matrix[x][y], matrix[x-1][y], matrix[x][y-1]);
		if(smaller != matrix[x][y]){
			if(matrix[x-1][y] != matrix[x][y-1]){
				if(smaller == matrix[x][y-1]){
					y -= 1;
				}else{
					x -= 1;
				}
			}
		}
	}

	if(error){
		res[0] = last_x;
		res[1] = last_y;
	}else{
		res[0] = --x;
		res[1] = --y;
	}
	return res;
}

static int* build_right_matrix(char *seq, char *motif, int mlen, int **matrix, Py_ssize_t start, int size, int max_error){
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

	for(x=1,y=1; x<=size && y<=size; x++,y++){
		h = seq[start+y];
		v = motif[(x-1)%mlen];
		
		//fill column, column number fixed
		if(i != y){
			for(i=1; i<x; i++){
				if(h == motif[(i-1)%mlen]){
					matrix[i][y] = matrix[i-1][y-1];
				}else{
					matrix[i][y] = min_three(matrix[i-1][y-1], matrix[i-1][y], matrix[i][y-1]) + 1;
				}
			}
		}
		//fill row, row number fixed
		if(j != x){
			for(j=1; j<y; j++){
				if(v == seq[start+j]){
					matrix[x][j] = matrix[x-1][j-1];
				}else{
					matrix[x][j] = min_three(matrix[x-1][j-1], matrix[x-1][j], matrix[x][j-1]) + 1;
				}
			}
		}

		i = y;
		j = x;

		if(h == v){
			matrix[x][y] = matrix[x-1][y-1];
			error = 0;
		}else{
			if(error == 0){
				last_x = x - 1;
				last_y = y - 1;
			}
			
			error++;

			if(error > max_error){
				break;
			}
			
			matrix[x][y] = min_three(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) + 1;
		}

		smaller = min_three(matrix[x][y], matrix[x-1][y], matrix[x][y-1]);
		if(smaller != matrix[x][y]){
			if(matrix[x-1][y] != matrix[x][y-1]){
				if(smaller == matrix[x][y-1]){
					y -= 1;
				}else{
					x -= 1;
				}
			}
		}
	}

	if(error){
		res[0] = last_x;
		res[1] = last_y;
	}else{
		res[0] = --x;
		res[1] = --y;
	}
	return res;
}

static int backtrace_matrix(int **matrix, int *diagonal, int *mat, int *sub, int *ins, int *del){
	int i = *diagonal;
	int j = *(diagonal+1);
	int cost;
	int r = j;

	while(i>0 && j>0){
		cost = min_three(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]);
		if(cost == matrix[i-1][j-1]){
			if(cost == matrix[i][j]){
				*mat += 1;
			}else{
				*sub += 1;
			}
			i--;
			j--;
		}else if(cost == matrix[i-1][j]){
			*del += 1;
			i--;
		}else{
			*ins += 1;
			j--;
		}
	}

	if(i>0){
		*del += 1;
	}else if(j>0){
		*ins += 1;
	}

	return r;
}

static PyObject* stripy_itrminer_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {"name", "seq", "min_motif_size", "max_motif_size", "seed_min_repeats", "seed_min_length", "max_continuous_errors", "max_extend_length", NULL};

	stripy_ITRMiner *obj = (stripy_ITRMiner *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	//initialize parameters
	obj->next_start = 0;
	obj->seed_minrep = 3;
	obj->seed_minlen = 8;
	obj->max_errors = 3;
	obj->max_motif = 6;
	obj->min_motif = 1;
	obj->max_extend_length = 2000;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiiiiii", keywords, &obj->seqname, &obj->seqobj, &obj->min_motif_size, &obj->max_motif_size, &obj->seed_minrep, &obj->seed_minlen, &obj->max_errors, &obj->max_extend_length)) {
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	return (PyObject *)self;
}

void stripy_itrminer_dealloc(stripy_ITRMiner *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* stripy_itrminer_repr(stripy_ITRMiner *self) {
	return PyUnicode_FromFormat("<ITRMiner> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* stripy_itrminer_as_list(stripy_ITRMiner *self) {
	PyObject *itrs = PyList_New(0);
	PyObject *tmp;
	
	Py_ssize_t seed_start;
	Py_ssize_t seed_end;
	int seed_length;
	int seed_repeat;
	int seed_good;

	int matches;
	int substitution;
	int insertion;
	int deletion;

	Py_ssize_t extend_start;
	int* extend_end;
	int extend_maxlen;
	int extend_len;
	int extend_total;

	Py_ssize_t tandem_start;
	Py_ssize_t tandem_end;
	int tandem_length;

	float align_rate;
	
	char* motif = (char *)malloc(self->max_motif + 1);

	int **matrix = initial_matrix(self->max_extend_length);

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

			if (see_repeat >= self->seed_min_repeats && seed_length >= self->seed_min_length) {
				const char *p = self->seq + seed_start;
				seed_good = 1;

				for (unsigned int k = 1; k < self->min_motif; ++k) {
					unsigned int l = 0;
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
					matches = seed_length;
					insertion = 0;
					deletion = 0;
					substitution = 0;
					extend_total = 0;

					//extend left flank
					extend_start = seed_start;
					extend_maxlen = extend_start;

					if (extend_maxlen > self->max_extend_length) {
						extend_maxlen = self->max_extend_length;
					}

					extend_end = build_left_matrix(self->seq, motif, j, matrix, extend_start, extend_maxlen, self->max_continuous_errors);
					extend_len = backtrace_matrix(matrix, extend_end, &matches, &substitution, &insertion, &deletion);
					extend_total += extend_len;
					
					tandem_start = extend_start - extend_len + 1;

					//extend right flank
					extend_start = seed_end;
					extend_maxlen = self->size - extend_start - 1;
					if (extend_maxlen > self->max_extend_length) {
						extend_maxlen = self->max_extend_length;
					}

					extend_end = build_right_matrix(self->self, motif, j, matrix, extend_start, extend_maxlen, self->max_continuous_errors);
					extend_len = backtrace_matrix(matrix, extend_end, &matches, &substitution, &insertion, &deletion);
					extend_total += extend_len;

					tandem_end = extend_start + extend_len + 1;
					tandem_length = tandem_end - tandem_start + 1;

					//calcuate the alignment rate of extended sequence
					if (extend_total > 0) {
						align_rate = (substitution*0.5 + insertion + deletion) / extend_total;
					} else {
						align_rate = 1;
					}

					if (align_rate >= 0.7) {
						tmp = Py_BuildValue("Onniiiii", self->seqname, tandem_start, tandem_end, tandem_length, matches, substitution, insertion, deletion);
						PyList_Append(itrs, tmp);
						Py_DECREF(tmp);

						i = tandem_end;
						break;
					}
				}
			}
			i = seed_start;
		}
	}
}

static PyMethodDef stripy_itrminer_methods[] = {
	{"as_list", (PyCFunction)stripy_itrminer_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject stripy_ITRMinerType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"ITRMiner",                        /* tp_name */
	sizeof(stripy_ITRMiner),          /* tp_basicsize */
	0,                              /* tp_itemsize */
	(destructor)stripy_itrminer_dealloc,   /* tp_dealloc */
	0,                              /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_reserved */
	(reprfunc)stripy_itrminer_repr,                              /* tp_repr */
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
	(getiterfunc)stripy_itrminer_iter,     /* tp_iter */
	(iternextfunc)stripy_itrminer_next,    /* tp_iternext */
	stripy_itrminer_methods,          /* tp_methods */
	0,          /* tp_members */
	0,                               /* tp_getset */
	0,                              /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	0,                              /* tp_init */
	0,            /* tp_alloc */
	stripy_itrminer_new,              /* tp_new */
};
