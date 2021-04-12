/*
 * tre.c -- Tandem Repeat Element
 * 
 * define tandem repeat element class object and methods
 *
 */

#define PY_SSIZE_T_CLEAN
#include "tre.h"
#include "structmember.h"

void stripy_tre_dealloc(stripy_TRE *self) {
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* stripy_tre_repr(stripy_TRE *self) {
	return PyUnicode_FromFormat("<%s> (%s)%d @ %s:%zd-%zd", Py_TYPE(self)->tp_name, self->motif, self->repeats, PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}

PyObject* stripy_tre_as_list(stripy_TRE *self) {
	return Py_BuildValue("OnnsII", self->seqid, self->start, self->end, self->motif, self->repeats, self->length);
}

PyObject* stripy_tre_as_string(stripy_TRE *self, PyObject *args, PyObject *kwargs) {
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
}

/* VNTR */
static PyMethodDef stripy_tre_methods[] = {
	{"as_list", (PyCFunction)stripy_tre_as_list, METH_NOARGS, NULL},
	{"as_string", (PyCFunction)stripy_tre_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyMemberDef stripy_tre_members[] = {
	{"seqid", T_OBJECT, offsetof(stripy_TRE, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(stripy_TRE, start), READONLY},
	{"end", T_PYSSIZET, offsetof(stripy_TRE, end), READONLY},
	{"motif", T_STRING, offsetof(stripy_TRE, motif), READONLY},
	{"repeats", T_UINT, offsetof(stripy_TRE, repeats), READONLY},
	{"length", T_UINT, offsetof(stripy_TRE, length), READONLY},
	{NULL}
};

PyTypeObject stripy_TREType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "TRE",                        /* tp_name */
    sizeof(stripy_TRE),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stripy_tre_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stripy_tre_repr,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                   /* tp_as_sequence */
    0,                   /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,             /* tp_flags */
    0,                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,     /* tp_iter */
    0,    /* tp_iternext */
    stripy_tre_methods,          /* tp_methods */
    stripy_tre_members,          /* tp_members */
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