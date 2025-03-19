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
/*
 * @param s str, motif sequence
 * @param l int, motif length
 * @param m int, min motif length
*/
static int is_redundant_motif(char *s, int l, int m) {
	int i, j, b;

	if (m == 1) {
		return 0;
	}

	for (j = 1; j < m; ++j) {
		if (l % j == 0) {
			b = l - j;
			i = 0;

			while ((i < b) && (s[i] == s[i+j])) {
				++i;
			}

			if (i == b) {
				return 1;
			}
		}
	}

	return 0;
}

static PyObject* pytrf_gtrfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	int i;

	static char* keywords[] = {"chrom", "seq", "min_motif", "max_motif", "min_repeat", "min_length", NULL};

	pytrf_GTRFinder *obj = (pytrf_GTRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	obj->min_motif = 1;
	obj->max_motif = 100;
	obj->min_repeat = 3;
	obj->min_length = 10;

	//initialize start search position
	obj->next_start = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiii", keywords, &obj->seqname, &obj->seqobj, &obj->min_motif, &obj->max_motif, &obj->min_repeat, &obj->min_length)) {
		return NULL;
	}

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);
	obj->motif = (char *)malloc(obj->max_motif + 1);

	obj->limit = (Py_ssize_t *)malloc(sizeof(Py_ssize_t) * (obj->max_motif+1));
	for (i = 0; i <= obj->max_motif; ++i) {
		obj->limit[i] = obj->size - i;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	return (PyObject *)obj;
}

static void pytrf_gtrfinder_dealloc(pytrf_GTRFinder *self) {
	if (self->limit) {
		free(self->limit);
	}

	free(self->motif);

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

	for (i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 'N' || self->seq[i] == 'n') {
			continue;
		}

		cs = i;
		for (j = self->min_motif; j <= self->max_motif; ++j) {
			while ((i < self->limit[j]) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			rl = i + j - cs;
			rn = rl / j;
			rl = rn * j;

			if (rn >= self->min_repeat && rl >= self->min_length) {
				memcpy(self->motif, self->seq+cs, j);
				self->motif[j] = '\0';

				if (is_redundant_motif(self->motif, j, self->min_motif)) {
					i = cs;
					continue;
				}

				pytrf_ETR *gtr = PyObject_New(pytrf_ETR, &pytrf_ETRType);

				gtr->mlen = j;
				gtr->repeat = rn;
				gtr->length = rl;
				gtr->start = cs + 1;
				gtr->end = cs + rl;
				gtr->seqid = Py_NewRef(self->seqname);
				gtr->seqobj = Py_NewRef(self->seqobj);
				gtr->motif = PyUnicode_FromString(self->motif);

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

	PyObject *gtrs = PyList_New(0);
	PyObject *tmp;

	for (i = 0; i < self->size; ++i) {
		if (self->seq[i] == 'N' || self->seq[i] == 'n') {
			continue;
		}

		cs = i;
		for (j = self->min_motif; j <= self->max_motif; ++j) {
			while ((i < self->limit[j]) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			rl = i + j - cs;
			rn = rl / j;
			rl = rn * j;

			if (rn >= self->min_repeat && rl >= self->min_length) {
				memcpy(self->motif, self->seq+cs, j);
				self->motif[j] = '\0';

				if (is_redundant_motif(self->motif, j, self->min_motif)) {
					i = cs;
					continue;
				}

				gs = cs + 1;
				ge = cs + rl;

				tmp = Py_BuildValue("Onnsiii", self->seqname, gs, ge, self->motif, j, rn, rl);
				PyList_Append(gtrs, tmp);
				Py_DECREF(tmp);
				
				i = ge - 1;
				break;
			}

			i = cs;
		}
	}

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
