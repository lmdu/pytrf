#define PY_SSIZE_T_CLEAN
#include "vntr.h"
#include "structmember.h"

static PyObject* stripy_vntrminer_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {"name", "seq", "min_motif", "max_motif", "min_repeats", NULL};

	stripy_VNTRMiner *obj = (stripy_VNTRMiner *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	obj->min_motif = 7;
	obj->max_motif = 30;
	obj->min_repeats = 2;

	//initialize start search position
	obj->next_start = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iii", keywords, &obj->seqname, &obj->seqobj, &obj->min_motif, &obj->max_motif, &obj->min_repeats)) {
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	return (PyObject *)obj;
}

void stripy_vntrminer_dealloc(stripy_VNTRMiner *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* stripy_vntrminer_repr(stripy_VNTRMiner *self) {
	return PyUnicode_FromFormat("<VNTRMiner> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* stripy_vntrminer_iter(stripy_VNTRMiner *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* stripy_vntrminer_next(stripy_VNTRMiner *self) {
	//current start position
	Py_ssize_t current_start;

	//repeat length
	unsigned int replen;

	//repeat number
	unsigned int repeats;

	//the motif is a right vntr motif
	int is_vntr;

	for (Py_ssize_t i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (unsigned int j = self->min_motif; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;
			repeats = replen/j;

			if (repeats >= self->min_repeats) {
				//check motif is real motif with length >= min motif size
				const char *p = self->seq+current_start;
				is_vntr = 1;
				for (unsigned int k = 1; k < self->min_motif; ++k) {
					int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l + k == j) {
						is_vntr = 0;
						break;
					}
				}

				//otherwise we found a vntr
				if (is_vntr) {
					stripy_VNTR *vntr = (stripy_VNTR *)PyObject_CallObject((PyObject *)&stripy_VNTRType, NULL);
					vntr->motif = (char *)malloc(j + 1);
					memcpy(vntr->motif, self->seq+current_start, j);
					vntr->motif[j] = '\0';
					vntr->seqid = self->seqname;
					Py_INCREF(vntr->seqid);
					vntr->repeats = repeats;
					vntr->length = repeats * j;
					vntr->start = current_start + 1;
					vntr->end = current_start + vntr->length;
					self->next_start = vntr->end;
					return (PyObject *)vntr;
				}
			}

			i = current_start;
		}
	}

	return NULL;
}

static PyObject* stripy_vntrminer_as_list(stripy_VNTRMiner *self) {
	PyObject *vntrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t current_start;
	unsigned int replen;
	unsigned int repeats;
	unsigned int length;
	int is_vntr;

	char* motif = (char *)malloc(self->max_motif + 1);

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (unsigned int j = self->min_motif; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;
			repeats = replen/j;

			if (repeats >= self->min_repeats) {
				//check motif is real motif with length >= min motif size
				const char *p = self->seq+current_start;
				is_vntr = 1;
				for (unsigned int k = 1; k < self->min_motif; ++k) {
					unsigned int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l + k == j) {
						is_vntr = 0;
						break;
					}
				}

				if (is_vntr) {
					memcpy(motif, self->seq+current_start, j);
					motif[j] = '\0';
					length = repeats * j;
					tmp = Py_BuildValue("OnnsII", self->seqname, current_start+1, current_start+length, motif, repeats, length);
					PyList_Append(vntrs, tmp);
					Py_DECREF(tmp);

					i -= replen % j;

					break;
				}
			}

			i = current_start;
		}
	}
	free(motif);
	return vntrs;
}

/* methods for SSR object */
/*void stripy_vntr_dealloc(stripy_VNTR *self) {
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* stripy_vntr_repr(stripy_VNTR *self) {
	return PyUnicode_FromFormat("<VNTR> (%s)%d @ %s:%zd-%zd", self->motif, self->repeats, PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}*/

static PyMethodDef stripy_vntrminer_methods[] = {
	{"as_list", (PyCFunction)stripy_vntrminer_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject stripy_VNTRMinerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "VNTRMiner",                        /* tp_name */
    sizeof(stripy_VNTRMiner),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stripy_vntrminer_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stripy_vntrminer_repr,                              /* tp_repr */
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
    (getiterfunc)stripy_vntrminer_iter,     /* tp_iter */
    (iternextfunc)stripy_vntrminer_next,    /* tp_iternext */
    stripy_vntrminer_methods,          /* tp_methods */
    0,          /* tp_members */
    0,                               /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    stripy_vntrminer_new,              /* tp_new */
};

/*PyObject* stripy_vntr_as_list(stripy_VNTR *self) {
	return Py_BuildValue("OnnsII", self->seqid, self->start, self->end, self->motif, self->repeats, self->length);
}

PyObject* stripy_vntr_as_string(stripy_VNTR *self, PyObject *args, PyObject *kwargs) {
	char *separator = "\t";
	int terminator = 0;
	static char* keywords[] = {"separator", "terminator", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|si", keywords, &separator, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%s%s%zd%s%zd%s%s%s%d%s%d%s", 
		PyUnicode_AsUTF8(self->seqid),
		separator,
		self->start,
		separator,
		self->end,
		separator,
		self->motif,
		separator,
		self->repeats,
		separator,
		self->length,
		terminator ? "\n" : ""
	);
}*/

/* VNTR */
/*static PyMethodDef stripy_vntr_methods[] = {
	{"as_list", (PyCFunction)stripy_vntr_as_list, METH_NOARGS, NULL},
	{"as_string", (PyCFunction)stripy_vntr_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyMemberDef stripy_vntr_members[] = {
	{"seqid", T_OBJECT, offsetof(stripy_VNTR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(stripy_VNTR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(stripy_VNTR, end), READONLY},
	{"motif", T_STRING, offsetof(stripy_VNTR, motif), READONLY},
	{"repeats", T_UINT, offsetof(stripy_VNTR, repeats), READONLY},
	{"length", T_UINT, offsetof(stripy_VNTR, length), READONLY},
	{NULL}
};*/

PyTypeObject stripy_VNTRType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "VNTR",                        /* tp_name */
    sizeof(stripy_VNTR),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    0,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    0,                              /* tp_repr */
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
    0,          /* tp_members */
    0,                               /* tp_getset */
    &stripy_TREType,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    PyType_GenericAlloc,            /* tp_alloc */
    PyType_GenericNew,              /* tp_new */
};
