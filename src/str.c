#define PY_SSIZE_T_CLEAN
#include "str.h"
#include "etr.h"
#include "structmember.h"

int pytrf_strfinder_set_min_repeats(pytrf_STRFinder *self, PyObject* minrep_obj) {
	if (minrep_obj) {
		if (PyList_Check(minrep_obj)) {
			minrep_obj = PyList_AsTuple(minrep_obj);
		}

		if (PyTuple_Check(minrep_obj)) {
			Py_ssize_t len = PyTuple_Size(minrep_obj);

			if (len != 6) {
				PyErr_SetString(PyExc_ValueError, "min_repeats list or tuple must contain six numbers");
				return 0;
			}

			for (Py_ssize_t i = 1; i < 7; ++i) {
				PyObject *p = PyTuple_GetItem(minrep_obj, i-1);
				if (PyLong_Check(p)) {
					self->min_lens[i] = PyLong_AsSsize_t(p) * i;
				} else {
					PyErr_SetString(PyExc_ValueError, "six number needed for min_repeats");
					return 0;
				}
			}
		} else if (PyDict_Check(minrep_obj)) {
			PyObject *key, *value;
			Py_ssize_t pos = 0;

			while (PyDict_Next(minrep_obj, &pos, &key, &value)) {
				if (PyLong_Check(value) || PyLong_Check(key)) {
					self->min_lens[PyLong_AsSsize_t(key)] = PyLong_AsSsize_t(value) * PyLong_AsSsize_t(key);
				} else {
					PyErr_SetString(PyExc_ValueError, "the key and value in min_repeats dict must be number");
					return 0;
				}
			}
		} else {
			PyErr_SetString(PyExc_TypeError, "min_repeats must be list, tuple or dict");
			return 0;
		}
	}

	return 1;
}

static PyObject* pytrf_strfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	PyObject *minrep_obj = NULL;
	static char* keywords[] = {"name", "seq", "min_repeats", NULL};

	pytrf_STRFinder *obj = (pytrf_STRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	//initialize start search position
	obj->next_start = 0;

	//initialize minimal repeats
	obj->min_lens[0] = 0;
	obj->min_lens[1] = 12 * 1;
	obj->min_lens[2] = 7 * 2;
	obj->min_lens[3] = 5 * 3;
	obj->min_lens[4] = 4 * 4;
	obj->min_lens[5] = 4 * 5;
	obj->min_lens[6] = 4 * 6;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|O", keywords, &obj->seqname, &obj->seqobj, &minrep_obj)) {
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	//parse minimal repeats
	if (!pytrf_strfinder_set_min_repeats(obj, minrep_obj)) {
		return NULL;
	}

	return (PyObject *)obj;
}

void pytrf_strfinder_dealloc(pytrf_STRFinder *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* pytrf_strfinder_repr(pytrf_STRFinder *self) {
	return PyUnicode_FromFormat("<STRFinder> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* pytrf_strfinder_iter(pytrf_STRFinder *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* pytrf_strfinder_next(pytrf_STRFinder *self) {
	//current slide position
	Py_ssize_t i;

	//current start position
	Py_ssize_t current_start;

	Py_ssize_t boundary;

	//motif length
	int j;

	//repeat length
	int replen;

	for (i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (j = 1; j <= 6; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= self->min_lens[j]) {
				pytrf_ETR *ssr = PyObject_New(pytrf_ETR, &pytrf_ETRType);
				ssr->motif = (char *)malloc(j + 1);
				memcpy(ssr->motif, self->seq+current_start, j);
				ssr->motif[j] = '\0';
				ssr->mlen = j;
				ssr->seqid = self->seqname;
				Py_INCREF(ssr->seqid);
				ssr->repeats = replen/j;
				ssr->length = ssr->repeats * j;
				ssr->start = current_start + 1;
				ssr->end = current_start + ssr->length;

				self->next_start = ssr->end;

				return (PyObject *)ssr;
			}
			
			i = current_start;
		}
	}

	return NULL;
}

static PyObject* pytrf_strfinder_reset_min_repeats(pytrf_STRFinder *self, PyObject *args, PyObject *kwargs) {
	PyObject *minrep_obj = NULL;
	static char* keywords[] = {"min_repeats", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", keywords, &minrep_obj)) {
		return NULL;
	}

	if (!pytrf_strfinder_set_min_repeats(self, minrep_obj)) {
		return NULL;
	}
	return NULL;
}

static PyObject* pytrf_strfinder_as_list(pytrf_STRFinder *self) {
	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t current_start;
	Py_ssize_t ssr_end;
	Py_ssize_t boundary;
	int replen;
	int repeats;
	int length;
	char *motif = (char *)malloc(7);

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = 1; j < 7; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= self->min_lens[j]) {
				memcpy(motif, self->seq+current_start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;
				ssr_end = current_start+length;
				tmp = Py_BuildValue("Onnsiii", self->seqname, current_start+1, ssr_end, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = ssr_end;
				break;
			} else {
				i = current_start;
			}
		}
	}

	free(motif);
	return ssrs;
}

static PyObject* pytrf_strfinder_as_test(pytrf_STRFinder *self) {
	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t current_start;
	Py_ssize_t ssr_end;
	Py_ssize_t boundary;
	int replen;
	int repeats;
	int length;
	//char *motif = (char *)malloc(7);
	char motif[7];

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = 1; j < 7; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= self->min_lens[j]) {
				memcpy(motif, self->seq+current_start, j);
				motif[j] = '\0';
				repeats = replen / j;
				length = repeats * j;
				ssr_end = current_start+length;
				tmp = Py_BuildValue("Onnsiii", self->seqname, current_start+1, ssr_end, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = ssr_end;
				break;
			} else {
				i = current_start;
			}
		}
	}

	//free(motif);
	return ssrs;
}

static PyMethodDef pytrf_strfinder_methods[] = {
	{"as_list", (PyCFunction)pytrf_strfinder_as_list, METH_NOARGS, NULL},
	{"as_test", (PyCFunction)pytrf_strfinder_as_test, METH_NOARGS, NULL},
	{"reset_min_repeats", (PyCFunction)pytrf_strfinder_reset_min_repeats, METH_VARARGS|METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject pytrf_STRFinderType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "STRFinder",
	.tp_basicsize = sizeof(pytrf_STRFinder),
	.tp_dealloc = (destructor)pytrf_strfinder_dealloc,
	.tp_repr = (reprfunc)pytrf_strfinder_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "short tandem repeat finder",
	.tp_iter = (getiterfunc)pytrf_strfinder_iter,
	.tp_iternext = (iternextfunc)pytrf_strfinder_next,
	.tp_methods = pytrf_strfinder_methods,
	.tp_new = pytrf_strfinder_new,
};
