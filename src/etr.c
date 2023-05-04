/*
 * etr.c -- Exact Tandem Repeat
 * 
 * define exact tandem repeat element class object and methods
 *
 */

#define PY_SSIZE_T_CLEAN
#include "etr.h"
#include "structmember.h"

void stria_etr_dealloc(stria_ETR *self) {
	free(self->motif);
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* stria_etr_repr(stria_ETR *self) {
	return PyUnicode_FromFormat("<ETR> (%s)%d @ %s:%zd-%zd", self->motif,
		self->repeats, PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}

PyObject* stria_etr_as_list(stria_ETR *self) {
	return Py_BuildValue("Onnsiii", self->seqid, self->start, self->end, self->motif,
		self->mlen, self->repeats, self->length);
}

PyObject* stria_etr_as_dict(stria_ETR *self) {
	return Py_BuildValue("{s:O,s:n,s:n,s:s,s:i,s:i,s:i}", "chrom", self->seqid,
						 "start", self->start, "end", self->end, "motif", self->motif,
						 "type", self->mlen, "repeats", self->repeats, "length", self->length);
}

PyObject* stria_etr_as_gff(stria_ETR *self) {
	return PyUnicode_FromFormat("%S\tstria\tETR\t%zd\t%zd\t.\t+\t.\tMotif=%s;Type=%d;Repeats=%d;Length=%d\n",
								self->seqid, self->start, self->end, self->motif, self->mlen, self->repeats,
								self->length);
}

PyObject* stria_etr_as_string(stria_ETR *self, PyObject *args, PyObject *kwargs) {
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

PyObject *stria_etr_get_seq(stria_ETR *self, void* closure) {
	PyObject* ret = PyUnicode_New(self->length, 127);
	Py_UCS1* p = PyUnicode_1BYTE_DATA(ret);
	memcpy(p, PyUnicode_AsUTF8(self->seqid)+self->start-1, self->length);
	return ret;
}

static PyMethodDef stria_etr_methods[] = {
	{"as_list", (PyCFunction)stria_etr_as_list, METH_NOARGS, NULL},
	{"as_dict", (PyCFunction)stria_etr_as_dict, METH_NOARGS, NULL},
	{"as_gff", (PyCFunction)stria_etr_as_gff, METH_NOARGS, NULL},
	{"as_string", (PyCFunction)stria_etr_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef stria_etr_getsets[] = {
	{"seq", (getter)stria_etr_get_seq, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef stria_etr_members[] = {
	{"chrom", T_OBJECT, offsetof(stria_ETR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(stria_ETR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(stria_ETR, end), READONLY},
	{"motif", T_STRING, offsetof(stria_ETR, motif), READONLY},
	{"type", T_INT, offsetof(stria_ETR, mlen), READONLY},
	{"repeats", T_INT, offsetof(stria_ETR, repeats), READONLY},
	{"length", T_INT, offsetof(stria_ETR, length), READONLY},
	{NULL}
};

PyTypeObject stria_ETRType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "ETR",
	.tp_basicsize = sizeof(stria_ETR),
	.tp_dealloc = (destructor)stria_etr_dealloc,
	.tp_repr = (reprfunc)stria_etr_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "tandem repeat element",
	.tp_methods = stria_etr_methods,
	.tp_members = stria_etr_members,
	.tp_getset = stria_etr_getsets,
	.tp_new = PyType_GenericNew,
};