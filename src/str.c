#define PY_SSIZE_T_CLEAN
#include "str.h"
#include "etr.h"
#include "compat.h"
#include "structmember.h"

static PyObject* pytrf_strfinder_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	int i;
	int mono = 12;
	int di = 7;
	int tri = 5;
	int tetra = 4;
	int penta = 4;
	int hexa = 4;

	pytrf_STRFinder *obj = (pytrf_STRFinder *)type->tp_alloc(type, 0);
	if (!obj) return NULL;
	
	static char* keywords[] = {"name", "seq", "mono", "di", "tri", "tetra", "penta", "hexa", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiiiii", keywords, &obj->seqname, &obj->seqobj, &mono, &di, &tri, &tetra, &penta, &hexa)) {
		return NULL;
	}

	//initialize start search position
	obj->next_start = 0;

	//initialize minimal repeats
	obj->min_lens[0] = 0;
	obj->min_lens[1] = mono * 1;
	obj->min_lens[2] = di * 2;
	obj->min_lens[3] = tri * 3;
	obj->min_lens[4] = tetra * 4;
	obj->min_lens[5] = penta * 5;
	obj->min_lens[6] = hexa * 6;

	obj->seq = PyUnicode_AsUTF8AndSize(obj->seqobj, &obj->size);

	for (i = 0; i < 7; ++i) {
		obj->boundary[i] = obj->size - i;
	}

	Py_INCREF(obj->seqname);
	Py_INCREF(obj->seqobj);

	return (PyObject *)obj;
}

static void pytrf_strfinder_dealloc(pytrf_STRFinder *self) {
	if (self->boundary) {
		free(self->boundary);
	}

	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* pytrf_strfinder_repr(pytrf_STRFinder *self) {
	return PyUnicode_FromFormat("<STRFinder> for sequence %S", self->seqname);
}

static PyObject* pytrf_strfinder_iter(pytrf_STRFinder *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* pytrf_strfinder_next(pytrf_STRFinder *self) {
	//current slide position
	Py_ssize_t i;

	//current start position
	Py_ssize_t cs;

	//boundary position
	Py_ssize_t b;

	//motif length
	int j;

	//repeat length
	int rl;

	for (i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		cs = i;
		for (j = 1; j <= 6; ++j) {
			b = self->boundary[j];

			while ((i < b) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			rl = i + j - cs;

			if (rl >= self->min_lens[j]) {
				pytrf_ETR *ssr = PyObject_New(pytrf_ETR, &pytrf_ETRType);

				ssr->mlen = j;
				ssr->repeat = rl/j;
				ssr->length = ssr->repeat * j;
				ssr->start = cs + 1;
				ssr->end = cs + ssr->length;
				ssr->seqid = Py_NewRef(self->seqname);
				ssr->seqobj = Py_NewRef(self->seqobj);
				ssr->motif = PyUnicode_Substring(self->seqobj, cs, cs + j);

				self->next_start = ssr->end;
				return (PyObject *)ssr;
			}

			i = cs;
		}
	}

	return NULL;
}

static PyObject* pytrf_strfinder_as_list(pytrf_STRFinder *self) {
	int j;

	//repeat length
	int rl;

	//ssr repeats
	int sr;

	//ssr length
	int sl;
	
	char motif[7];

	Py_ssize_t i;
	
	//current start position
	Py_ssize_t cs;
	
	//str end position
	Py_ssize_t se;

	//boundary
	Py_ssize_t b;

	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;

	for (i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		cs = i;
		for (j = 1; j < 7; ++j) {
			b = self->boundary[j];

			while ((i < b) && (self->seq[i] == self->seq[i+j])) {
				++i;
			}

			rl = i + j - cs;

			if (rl >= self->min_lens[j]) {
				memcpy(motif, self->seq+cs, j);
				motif[j] = '\0';
				sr = rl / j;
				sl = sr * j;
				se = cs + sl;
				tmp = Py_BuildValue("Onnsiii", self->seqname, cs+1, se, motif, j, sr, sl);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = se - 1;
				break;
			}
			i = cs;
		}
	}

	return ssrs;
}

static PyMethodDef pytrf_strfinder_methods[] = {
	{"as_list", (PyCFunction)pytrf_strfinder_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject pytrf_STRFinderType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "STRFinder",
	.tp_basicsize = sizeof(pytrf_STRFinder),
	.tp_dealloc = (destructor)pytrf_strfinder_dealloc,
	.tp_repr = (reprfunc)pytrf_strfinder_repr,
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_doc = "short tandem repeat finder",
	.tp_iter = (getiterfunc)pytrf_strfinder_iter,
	.tp_iternext = (iternextfunc)pytrf_strfinder_next,
	.tp_methods = pytrf_strfinder_methods,
	.tp_new = pytrf_strfinder_new,
};
