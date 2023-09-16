/*
 * etr.c -- Exact Tandem Repeat
 * 
 * define exact tandem repeat element class object and methods
 *
 */

#define PY_SSIZE_T_CLEAN
#include "etr.h"
#include "structmember.h"

static void pytrf_etr_dealloc(pytrf_ETR *self) {
	Py_DECREF(self->motif);
	Py_DECREF(self->seqid);
	Py_DECREF(self->seqobj);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* pytrf_etr_repr(pytrf_ETR *self) {
	return PyUnicode_FromFormat("<ETR> (%S)%d @ %S:%zd-%zd", self->motif,
		self->repeats, self->seqid, self->start, self->end);
}

static PyObject* pytrf_etr_as_list(pytrf_ETR *self) {
	return Py_BuildValue("OnnOiii", self->seqid, self->start, self->end, self->motif,
		self->mlen, self->repeats, self->length);
}

static PyObject* pytrf_etr_as_dict(pytrf_ETR *self) {
	return Py_BuildValue("{s:O,s:n,s:n,s:O,s:i,s:i,s:i}", "chrom", self->seqid,
						 "start", self->start, "end", self->end, "motif", self->motif,
						 "type", self->mlen, "repeats", self->repeats, "length", self->length);
}

static PyObject* pytrf_etr_as_gff(pytrf_ETR *self, PyObject *args, PyObject *kwargs) {
	char *terminator = "";
	static char* keywords[] = {"terminator", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|s", keywords, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%S\tpytrf\tETR\t%zd\t%zd\t.\t+\t.\tMotif=%S;Type=%d;Repeats=%d;Length=%d%s",
								self->seqid, self->start, self->end, self->motif, self->mlen, self->repeats,
								self->length, terminator);
}

static PyObject* pytrf_etr_as_string(pytrf_ETR *self, PyObject *args, PyObject *kwargs) {
	char *separator = "\t";
	char *terminator = "";
	static char* keywords[] = {"separator", "terminator", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|ss", keywords, &separator, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%S%s%zd%s%zd%s%S%s%d%s%d%s%d%s", self->seqid, separator,
								self->start, separator, self->end, separator, self->motif,
								separator, self->mlen, separator, self->repeats, separator,
								self->length, terminator);
}

static PyObject* pytrf_etr_get_seq(pytrf_ETR *self, void* closure) {
	return PyUnicode_Substring(self->seqobj, self->start-1, self->end);
}

static PyMethodDef pytrf_etr_methods[] = {
	{"as_list", (PyCFunction)pytrf_etr_as_list, METH_NOARGS, NULL},
	{"as_dict", (PyCFunction)pytrf_etr_as_dict, METH_NOARGS, NULL},
	{"as_gff", (PyCFunction)pytrf_etr_as_gff, METH_VARARGS | METH_KEYWORDS, NULL},
	{"as_string", (PyCFunction)pytrf_etr_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef pytrf_etr_getsets[] = {
	{"seq", (getter)pytrf_etr_get_seq, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pytrf_etr_members[] = {
	{"chrom", T_OBJECT, offsetof(pytrf_ETR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(pytrf_ETR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(pytrf_ETR, end), READONLY},
	{"motif", T_OBJECT, offsetof(pytrf_ETR, motif), READONLY},
	{"type", T_INT, offsetof(pytrf_ETR, mlen), READONLY},
	{"repeats", T_INT, offsetof(pytrf_ETR, repeats), READONLY},
	{"length", T_INT, offsetof(pytrf_ETR, length), READONLY},
	{NULL}
};

PyTypeObject pytrf_ETRType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "ETR",
	.tp_basicsize = sizeof(pytrf_ETR),
	.tp_dealloc = (destructor)pytrf_etr_dealloc,
	.tp_repr = (reprfunc)pytrf_etr_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "exact or perfect tandem repeat",
	.tp_methods = pytrf_etr_methods,
	.tp_members = pytrf_etr_members,
	.tp_getset = pytrf_etr_getsets,
	.tp_new = PyType_GenericNew,
};