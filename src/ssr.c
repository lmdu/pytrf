#define PY_SSIZE_T_CLEAN
#include "ssr.h"
#include "etr.h"
#include "structmember.h"

int stria_ssrminer_set_min_repeats(stria_SSRMiner *self, PyObject* minrep_obj) {
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

static PyObject* stria_ssrminer_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	PyObject *minrep_obj = NULL;
	static char* keywords[] = {"name", "seq", "min_repeats", NULL};

	stria_SSRMiner *obj = (stria_SSRMiner *)type->tp_alloc(type, 0);
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
	if (!stria_ssrminer_set_min_repeats(obj, minrep_obj)) {
		return NULL;
	}

	return (PyObject *)obj;
}

void stria_ssrminer_dealloc(stria_SSRMiner *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* stria_ssrminer_repr(stria_SSRMiner *self) {
	return PyUnicode_FromFormat("<SSRMiner> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* stria_ssrminer_iter(stria_SSRMiner *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* stria_ssrminer_next(stria_SSRMiner *self) {
	//current slide position
	Py_ssize_t i;

	//current start position
	Py_ssize_t current_start;

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
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= self->min_lens[j]) {
				stria_ETR *ssr = PyObject_New(stria_ETR, &stria_ETRType);
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

static PyObject* stria_ssrminer_reset_min_repeats(stria_SSRMiner *self, PyObject *args, PyObject *kwargs) {
	PyObject *minrep_obj = NULL;
	static char* keywords[] = {"min_repeats", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", keywords, &minrep_obj)) {
		return NULL;
	}

	if (!stria_ssrminer_set_min_repeats(self, minrep_obj)) {
		return NULL;
	}
	return NULL;
}

static PyObject* stria_ssrminer_as_list(stria_SSRMiner *self) {
	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t current_start;
	int replen;
	int repeats;
	int length;
	char motif[7];

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = 1; j < 7; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= self->min_lens[j]) {
				memcpy(motif, self->seq+current_start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;

				tmp = Py_BuildValue("Onnsiii", self->seqname, current_start+1, current_start+length, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i -= replen % j;

				break;
			} else {
				i = current_start;
			}
		}
	}

	return ssrs;
}

static PyMethodDef stria_ssrminer_methods[] = {
	{"as_list", (PyCFunction)stria_ssrminer_as_list, METH_NOARGS, NULL},
	{"reset_min_repeats", (PyCFunction)stria_ssrminer_reset_min_repeats, METH_VARARGS|METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject stria_SSRMinerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "SSRMiner",                        /* tp_name */
    sizeof(stria_SSRMiner),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stria_ssrminer_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stria_ssrminer_repr,                              /* tp_repr */
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
    (getiterfunc)stria_ssrminer_iter,     /* tp_iter */
    (iternextfunc)stria_ssrminer_next,    /* tp_iternext */
    stria_ssrminer_methods,          /* tp_methods */
    0,          /* tp_members */
    0,                               /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    0,            /* tp_alloc */
    stria_ssrminer_new,              /* tp_new */
};
