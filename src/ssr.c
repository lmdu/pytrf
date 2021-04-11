#define PY_SSIZE_T_CLEAN
#include "ssr.h"
#include "structmember.h"

int stripy_ssrminer_set_min_repeats(stripy_SSRMiner *self, PyObject* minrep_obj) {
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

static PyObject* stripy_ssrminer_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	PyObject *minrep_obj = NULL;
	static char* keywords[] = {"name", "seq", "min_repeats", NULL};

	stripy_SSRMiner *obj = (stripy_SSRMiner *)type->tp_alloc(type, 0);
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
	if (!stripy_ssrminer_set_min_repeats(obj, minrep_obj)) {
		return NULL;
	}

	return (PyObject *)obj;
}

void stripy_ssrminer_dealloc(stripy_SSRMiner *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* stripy_ssrminer_repr(stripy_SSRMiner *self) {
	return PyUnicode_FromFormat("<SSRMiner> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* stripy_ssrminer_iter(stripy_SSRMiner *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* stripy_ssrminer_next(stripy_SSRMiner *self) {
	//current slide position
	Py_ssize_t i;

	//current start position
	Py_ssize_t current_start;

	//motif length
	unsigned int j;

	//repeat length
	unsigned int replen;

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
				stripy_SSR *ssr = (stripy_SSR *)PyObject_CallObject((PyObject *)&stripy_SSRType, NULL);
				memcpy(ssr->motif, self->seq+self->next_start, j);
				ssr->motif[j] = '\0';
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

static PyObject* stripy_ssrminer_reset_min_repeats(stripy_SSRMiner *self, PyObject *args, PyObject *kwargs) {
	PyObject *minrep_obj = NULL;
	static char* keywords[] = {"min_repeats", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", keywords, &minrep_obj)) {
		return NULL;
	}

	if (!stripy_ssrminer_set_min_repeats(self, minrep_obj)) {
		return NULL;
	}
	return NULL;
}

static PyObject* stripy_ssrminer_as_list(stripy_SSRMiner *self) {
	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t start;
	unsigned int replen;
	unsigned int repeats;
	unsigned int length;
	char motif[7];

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		start = i;
		for (int j = 1; j < 7; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - start;

			if (replen >= self->min_lens[j]) {
				memcpy(motif, self->seq+start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;

				tmp = Py_BuildValue("OnnsII", self->seqname, start+1, start+length, motif, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i -= replen % j;

				break;
			} else {
				i = start;
			}
		}
	}

	return ssrs;
}

/* methods for SSR object */
void stripy_ssr_dealloc(stripy_SSR *self) {
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* stripy_ssr_repr(stripy_SSR *self) {
	return PyUnicode_FromFormat("<SSR> (%s)%d @ %s:%zd-%zd", self->motif, self->repeats, PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}

static PyMethodDef stripy_ssrminer_methods[] = {
	{"as_list", (PyCFunction)stripy_ssrminer_as_list, METH_NOARGS, NULL},
	{"reset_min_repeats", (PyCFunction)stripy_ssrminer_reset_min_repeats, METH_VARARGS|METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject stripy_SSRMinerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "SSRMiner",                        /* tp_name */
    sizeof(stripy_SSRMiner),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stripy_ssrminer_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stripy_ssrminer_repr,                              /* tp_repr */
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
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)stripy_ssrminer_iter,     /* tp_iter */
    (iternextfunc)stripy_ssrminer_next,    /* tp_iternext */
    stripy_ssrminer_methods,          /* tp_methods */
    0,          /* tp_members */
    0,                               /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    stripy_ssrminer_new,              /* tp_new */
};

static PyMemberDef stripy_ssr_members[] = {
	{"seqid", T_OBJECT, offsetof(stripy_SSR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(stripy_SSR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(stripy_SSR, end), READONLY},
	{"motif", T_STRING, offsetof(stripy_SSR, motif), READONLY},
	{"repeats", T_UINT, offsetof(stripy_SSR, repeats), READONLY},
	{"length", T_UINT, offsetof(stripy_SSR, length), READONLY},
	{NULL}
};

PyTypeObject stripy_SSRType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "SSR",                        /* tp_name */
    sizeof(stripy_SSR),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stripy_ssr_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stripy_ssr_repr,                              /* tp_repr */
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
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,     /* tp_iter */
    0,    /* tp_iternext */
    0,          /* tp_methods */
    stripy_ssr_members,          /* tp_members */
    0,                               /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    PyType_GenericNew,              /* tp_new */
};
