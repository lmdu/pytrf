/* 
 * gtr.c -- Generic tandem repeat finder
 *
 * define the finder class for generic tandem repeat
 *
 */

#define PY_SSIZE_T_CLEAN
#include "gtr.h"
#include "etr.h"
#include "structmember.h"

PyObject* pytrf_gtrfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {"chrom", "seq", "min_motif_size", "max_motif_size", "min_repeat", NULL};

	pytrf_GTRFinder *obj = (pytrf_GTRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	obj->min_motif = 10;
	obj->max_motif = 100;
	obj->min_repeat = 3;

	//initialize start search position
	obj->next_start = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iii", keywords, &obj->seqname, &obj->seqobj, &obj->min_motif, &obj->max_motif, &obj->min_repeat)) {
		return NULL;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	return (PyObject *)obj;
}

void pytrf_gtrfinder_dealloc(pytrf_GTRFinder *self) {
	self->seq = NULL;
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* pytrf_gtrfinder_repr(pytrf_GTRFinder *self) {
	return PyUnicode_FromFormat("<GTRFinder> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* pytrf_gtrfinder_iter(pytrf_GTRFinder *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* pytrf_gtrfinder_next(pytrf_GTRFinder *self) {
	//current start position
	Py_ssize_t current_start;

	//boundary
	Py_ssize_t boundary;

	//repeat length
	int replen;

	//repeat number
	int repeats;

	//the motif is a right gtr motif
	int is_gtr;

	for (Py_ssize_t i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;
			repeats = replen/j;

			if (repeats >= self->min_repeat) {
				//check motif is real motif with length >= min motif size
				const char *p = self->seq+current_start;
				is_gtr = 1;
				for (int k = 1; k < self->min_motif; ++k) {
					int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l + k == j) {
						is_gtr = 0;
						break;
					}
				}

				//otherwise we found a gtr
				if (is_gtr) {
					//stria_VNTR *gtr = (stria_VNTR *)PyObject_CallObject((PyObject *)&stria_VNTRType, NULL);
					pytrf_ETR *gtr = PyObject_New(pytrf_ETR, &pytrf_ETRType);
					gtr->motif = (char *)malloc(j + 1);
					memcpy(gtr->motif, self->seq+current_start, j);
					gtr->motif[j] = '\0';
					gtr->mlen = j;
					gtr->seqid = self->seqname;
					Py_INCREF(gtr->seqid);
					gtr->repeats = repeats;
					gtr->length = repeats * j;
					gtr->start = current_start + 1;
					gtr->end = current_start + gtr->length;
					self->next_start = gtr->end;
					return (PyObject *)gtr;
				}
			}

			i = current_start;
		}
	}

	return NULL;
}

static PyObject* pytrf_gtrfinder_as_list(pytrf_GTRFinder *self) {
	PyObject *gtrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t current_start;
	Py_ssize_t gtr_end;
	Py_ssize_t boundary;
	int repeats;
	int length;

	char* motif = (char *)malloc(self->max_motif + 1);

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			boundary = self->size - j;

			while ((i < boundary) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			length = i + j - current_start;
			repeats = length/j;

			if (repeats >= self->min_repeat) {
				//check motif is real motif with length >= min motif size
				const char *p = self->seq+current_start;
				int is_gtr = 1;
				int l = 0;

				for (int k = 1; k < self->min_motif; ++k) {
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l + k == j) {
						is_gtr = 0;
						break;
					}
				}

				if (is_gtr) {
					memcpy(motif, self->seq+current_start, j);
					motif[j] = '\0';
					length = repeats * j;
					gtr_end = current_start+length;
					tmp = Py_BuildValue("Onnsiii", self->seqname, current_start+1, gtr_end,
										motif, j, repeats, length);
					PyList_Append(gtrs, tmp);
					Py_DECREF(tmp);

					i = gtr_end;
					break;
				}
			}

			i = current_start;
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
