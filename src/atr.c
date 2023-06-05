/*
 * define class for approximate tandem repeat element
 */
#define PY_SSIZE_T_CLEAN
#include "atr.h"
#include "structmember.h"

void pytrf_atr_dealloc(pytrf_ATR *self) {
	free(self->motif);
	Py_DECREF(self->seqid);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

PyObject* pytrf_atr_repr(pytrf_ATR *self) {
	return PyUnicode_FromFormat("<ATR> %s @ %s:%zd-%zd", self->motif,
								PyUnicode_AsUTF8(self->seqid), self->start, self->end);
}

PyObject* pytrf_atr_as_list(pytrf_ATR *self) {
	return Py_BuildValue("Onnsiiiiiif", self->seqid, self->start, self->end, self->motif, self->mlen,
						self->length, self->matches, self->substitutions, self->insertions,
						self->deletions, self->identity);
}

PyObject* pytrf_atr_as_dict(pytrf_ATR *self) {
	return Py_BuildValue("{s:O,s:n,s:n,s:s,s:i,s:i,s:i,s:i,s:i,s:i,s:d}", "chrom", self->seqid,
						 "start", self->start, "end", self->end, "motif", self->motif, "type", self->mlen,
						 "length", self->length, "matches", self->matches, "substitutions", self->substitutions,
						 "insertions", self->insertions, "deletions", self->deletions, "identity", self->identity);
}

PyObject* pytrf_atr_as_gff(pytrf_ATR *self, PyObject *args, PyObject *kwargs) {
	char *terminator = "";
	static char* keywords[] = {"terminator", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|s", keywords, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%S\tpytrf\tATR\t%zd\t%zd\t.\t+\t.\tMotif=%s;Type=%d;Length=%d;Match=%d;"
								"Substitutions=%d;Insertions=%d;Deletions=%d;Identity=%R%s", self->seqid, self->start,
								self->end, self->motif, self->mlen, self->length, self->matches, self->substitutions, 
								self->insertions, self->deletions, PyFloat_FromDouble(self->identity), terminator);
}

PyObject* pytrf_atr_as_string(pytrf_ATR *self, PyObject *args, PyObject *kwargs) {
	char *separator = "\t";
	char *terminator = "";

	static char* keywords[] = {"separator", "terminator", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|ss", keywords, &separator, &terminator)) {
		return NULL;
	}

	return PyUnicode_FromFormat("%S%s%zd%s%zd%s%s%s%d%s%d%s%d%s%d%s%d%s%d%s%R%s", self->seqid, separator,
								self->start, separator, self->end, separator, self->motif, separator,
								self->mlen, separator, self->length, separator, self->matches, separator,
								self->substitutions, separator, self->insertions, separator, self->deletions,
								separator, PyFloat_FromDouble(self->identity), terminator);
}

PyObject *pytrf_atr_get_seq(pytrf_ATR *self, void* closure) {
	PyObject* ret = PyUnicode_New(self->length, 127);
	Py_UCS1* p = PyUnicode_1BYTE_DATA(ret);
	memcpy(p, PyUnicode_AsUTF8(self->seqid)+self->start-1, self->length);
	return ret;
}

static PyMethodDef pytrf_atr_methods[] = {
	{"as_list", (PyCFunction)pytrf_atr_as_list, METH_NOARGS, NULL},
	{"as_dict", (PyCFunction)pytrf_atr_as_dict, METH_NOARGS, NULL},
	{"as_gff", (PyCFunction)pytrf_atr_as_gff, METH_NOARGS, NULL},
	{"as_string", (PyCFunction)pytrf_atr_as_string, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static PyGetSetDef pytrf_atr_getsets[] = {
	{"seq", (getter)pytrf_atr_get_seq, NULL, NULL, NULL},
	{NULL}
};

static PyMemberDef pytrf_atr_members[] = {
	{"chrom", T_OBJECT, offsetof(pytrf_ATR, seqid), READONLY},
	{"start", T_PYSSIZET, offsetof(pytrf_ATR, start), READONLY},
	{"end", T_PYSSIZET, offsetof(pytrf_ATR, end), READONLY},
	{"motif", T_STRING, offsetof(pytrf_ATR, motif), READONLY},
	{"type", T_INT, offsetof(pytrf_ATR, mlen), READONLY},
	{"length", T_INT, offsetof(pytrf_ATR, length), READONLY},
	{"matches", T_INT, offsetof(pytrf_ATR, matches), READONLY},
	{"substitutions", T_INT, offsetof(pytrf_ATR, substitutions), READONLY},
	{"insertions", T_INT, offsetof(pytrf_ATR, insertions), READONLY},
	{"deletions", T_INT, offsetof(pytrf_ATR, deletions), READONLY},
	{"identity", T_DOUBLE, offsetof(pytrf_ATR, identity), READONLY},
	{NULL}
};

PyTypeObject pytrf_ATRType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "ATR",
	.tp_basicsize = sizeof(pytrf_ATR),
	.tp_dealloc = (destructor)pytrf_atr_dealloc,
	.tp_repr = (reprfunc)pytrf_atr_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "approximate tandem repeat element",
	.tp_methods = pytrf_atr_methods,
	.tp_members = pytrf_atr_members,
	.tp_getset = pytrf_atr_getsets,
	.tp_new = PyType_GenericNew,
};