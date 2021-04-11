#define PY_SSIZE_T_CLEAN
#include "ssr.h"
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
	//current slide position
	Py_ssize_t i;

	//current start position
	Py_ssize_t current_start;

	//motif length
	unsigned int j;

	//repeat length
	unsigned int replen;

	//repeat number
	unsigned int repeats;

	for (i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (j = 1; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;


			if (replen >= self->min_lens[j]) {
				stripy_VNTR *vntr = (stripy_VNTR *)PyObject_CallObject((PyObject *)&stripy_VNTRType, NULL);
				vntr->motif = (char *)malloc(j + 1);
				memcpy(vntr->motif, self->seq+self->next_start, j);
				vntr->motif[j] = '\0';
				vntr->seqid = self->seqname;
				Py_INCREF(vntr->seqid);
				vntr->repeats = replen/j;
				vntr->length = vntr->repeats * j;
				vntr->start = current_start + 1;
				vntr->end = current_start + vntr->length;

				self->next_start = vntr->end;

				return (PyObject *)vntr;
			}
			
			i = current_start;
		}
	}

	return NULL;
}

static PyObject* stripy_vntrminer_as_list(stripy_VNTRMiner *self) {
	PyObject *vntrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t start;
	unsigned int replen;
	unsigned int repeats;
	unsigned int length;
	char* motif = (char *)malloc(self->max_motif + 1);

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
