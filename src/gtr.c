/* 
 * gtr.c -- Generic tandem repeat finder
 *
 * define the finder class for generic tandem repeat
 *
 */

#define PY_SSIZE_T_CLEAN
#include "gtr.h"
#include "etr.h"
#include "compat.h"
#include "structmember.h"

static PyObject* pytrf_gtrfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {"chrom", "seq", "max_motif", "min_repeat", "min_length", NULL};

	pytrf_GTRFinder *obj = (pytrf_GTRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	obj->max_motif = 30;
	obj->min_repeat = 3;
	obj->min_length = 10;

	//initialize start search position
	obj->next_start = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iii", keywords, &obj->seqname, &obj->seqobj, &obj->max_motif, &obj->min_repeat, &obj->min_length)) {
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	return (PyObject *)obj;
}

static void pytrf_gtrfinder_dealloc(pytrf_GTRFinder *self) {
	self->seq = NULL;
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* pytrf_gtrfinder_repr(pytrf_GTRFinder *self) {
	return PyUnicode_FromFormat("<GTRFinder> for sequence %S", self->seqname);
}

static PyObject* pytrf_gtrfinder_iter(pytrf_GTRFinder *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* pytrf_gtrfinder_next(pytrf_GTRFinder *self) {
	int j;
	Py_ssize_t i;

	//repeat length
	int rl;

	//repeat number
	int rn;

	//current start position
	Py_ssize_t cs;

	//end position
	Py_ssize_t ep;

	for (i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		cs = i;
		for (j = 1; j <= self->max_motif; ++j) {
			ep = self->size - j;

			while ((i < ep) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			rl = i + j - cs;
			rn = rl/j;
			rl = rn*j;

			if (rn >= self->min_repeat && rl >= self->min_length) {
				pytrf_ETR *gtr = PyObject_New(pytrf_ETR, &pytrf_ETRType);
				gtr->mlen = j;
				gtr->repeats = rn;
				gtr->length = rl;
				gtr->start = cs + 1;
				gtr->end = cs + rl;
				gtr->seqid = Py_NewRef(self->seqname);
				gtr->seqobj = Py_NewRef(self->seqobj);
				gtr->motif = PyUnicode_Substring(self->seqobj, cs, cs + j);

				self->next_start = gtr->end;
				return (PyObject *)gtr;
			}

			i = cs;
		}
	}

	return NULL;
}

static PyObject* pytrf_gtrfinder_as_list(pytrf_GTRFinder *self) {
	int j;
	Py_ssize_t i;

	//repeat number
	int rn;

	//repeat length
	int rl;

	//current start position
	Py_ssize_t cs;

	//gtr start and end position
	Py_ssize_t gs;
	Py_ssize_t ge;

	//end position
	Py_ssize_t ep;

	//motif cache
	char *motif;

	PyObject *gtrs = PyList_New(0);
	PyObject *tmp;

	motif = (char *)malloc(self->max_motif + 1);

	for (i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		cs = i;
		for (j = 1; j <= self->max_motif; ++j) {
			ep = self->size - j;

			while ((i < ep) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			rl = i + j - cs;
			rn = rl / j;
			rl = rn * j;

			if (rn >= self->min_repeat && rl >= self->min_length) {
				memcpy(motif, self->seq+cs, j);
				motif[j] = '\0';
				gs = cs + 1;
				ge = cs + rl;
				tmp = Py_BuildValue("Onnsiii", self->seqname, gs, ge, motif, j, rn, rl);
				PyList_Append(gtrs, tmp);
				Py_DECREF(tmp);
				i = ge;
				break;
			}

			i = cs;
		}
	}

	free(motif);
	return gtrs;
}

static PyMethodDef pytrf_gtrfinder_methods[] = {
	{"as_list", (PyCFunction)pytrf_gtrfinder_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject pytrf_GTRFinderType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "GTRFinder",
    .tp_basicsize = sizeof(pytrf_GTRFinder),
    .tp_dealloc = (destructor)pytrf_gtrfinder_dealloc,
    .tp_repr = (reprfunc)pytrf_gtrfinder_repr,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "generic tandem repeat finder",
    .tp_iter = (getiterfunc)pytrf_gtrfinder_iter,
    .tp_iternext = (iternextfunc)pytrf_gtrfinder_next,
    .tp_methods = pytrf_gtrfinder_methods,
    .tp_new = pytrf_gtrfinder_new,
};
