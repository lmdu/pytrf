/*
 * etr.c -- Exact Tandem Repeat
 * 
 * define exact tandem repeat element class object and methods
 *
 */

#define PY_SSIZE_T_CLEAN
#include "etr.h"
#include "structmember.h"

void stripy_etr_dealloc(stripy_ETR *self) {
	free(self->motif);
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* stripy_etr_repr(stripy_ETR *self) {
	return PyUnicode_FromFormat("<ETR> (%s)%d @ %s:%zd-%zd", self->motif,
		self->repeats, PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}

PyObject* stripy_etr_as_list(stripy_ETR *self) {
	return Py_BuildValue("Onnsiii", self->seqid, self->start, self->end, self->motif,
		self->mlen, self->repeats, self->length);
}

PyObject* stripy_etr_as_string(stripy_ETR *self, PyObject *args, PyObject *kwargs) {
	char *separator = "\t";
	char *terminator = "";
	static char* keywords[] = {"separator", "terminator", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|ss", keywords, &separator, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%S%s%zd%s%zd%s%s%s%d%s%d%s%d%s", self->seqid, separator,
								self->start, separator, self->end, separator, self->motif,
								separator, self->mlen, separator, self->repeats, separator,
								self->length, terminator);
}

PyObject *stripy_etr_get_seq(stripy_ETR *self, void* closure) {
	PyObject* ret = PyUnicode_New(self->length, 127);
	Py_UCS1* p = PyUnicode_1BYTE_DATA(ret);
	memcpy(p, PyUnicode_AsUTF8(self->seqid)+self->start-1, self->length);
	return ret;
}

static PyMethodDef stripy_etr_methods[] = {
	{"as_list", (PyCFunction)stripy_etr_as_list, METH_NOARGS, NULL},
	{"as_string", (PyCFunction)stripy_etr_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef stripy_etr_getsets[] = {
	{"seq", (getter)stripy_etr_get_seq, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef stripy_etr_members[] = {
	{"seqid", T_OBJECT, offsetof(stripy_ETR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(stripy_ETR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(stripy_ETR, end), READONLY},
	{"motif", T_STRING, offsetof(stripy_ETR, motif), READONLY},
	{"type", T_INT, offsetof(stripy_ETR, mlen), READONLY},
	{"repeats", T_INT, offsetof(stripy_ETR, repeats), READONLY},
	{"length", T_INT, offsetof(stripy_ETR, length), READONLY},
	{NULL}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
};

PyTypeObject stripy_ETRType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "ETR",                        /* tp_name */
    sizeof(stripy_ETR),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stripy_etr_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stripy_etr_repr,                              /* tp_repr */
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
    stripy_etr_methods,          /* tp_methods */
    stripy_etr_members,          /* tp_members */
    stripy_etr_getsets,                               /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    0,            /* tp_alloc */
    PyType_GenericNew,              /* tp_new */
};