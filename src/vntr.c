#define PY_SSIZE_T_CLEAN
#include "vntr.h"
#include "etr.h"
#include "structmember.h"

static PyObject* stria_vntrminer_new(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
	static char* keywords[] = {"name", "seq", "min_motif_size", "max_motif_size", "min_repeat", NULL};

	stria_VNTRMiner *obj = (stria_VNTRMiner *)type->tp_alloc(type, 0);
	if (!obj) return NULL;

	obj->min_motif = 7;
	obj->max_motif = 30;
	obj->min_repeat = 2;

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

void stria_vntrminer_dealloc(stria_VNTRMiner *self) {
	Py_DECREF(self->seqname);
	Py_DECREF(self->seqobj);
	self->seq = NULL;
	Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject* stria_vntrminer_repr(stria_VNTRMiner *self) {
	return PyUnicode_FromFormat("<VNTRMiner> for sequence %s", PyUnicode_AsUTF8(self->seqname));
}

static PyObject* stria_vntrminer_iter(stria_VNTRMiner *self) {
	self->next_start = 0;
	Py_INCREF(self);
	return (PyObject *)self;
}

static PyObject* stria_vntrminer_next(stria_VNTRMiner *self) {
	//current start position
	Py_ssize_t current_start;

	//repeat length
	int replen;

	//repeat number
	int repeats;

	//the motif is a right vntr motif
	int is_vntr;

	for (Py_ssize_t i = self->next_start; i < self->size; ++i) {
		//remove unkown base
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;
			repeats = replen/j;

			if (repeats >= self->min_repeat) {
				//check motif is real motif with length >= min motif size
				const char *p = self->seq+current_start;
				is_vntr = 1;
				for (int k = 1; k < self->min_motif; ++k) {
					int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l + k == j) {
						is_vntr = 0;
						break;
					}
				}

				//otherwise we found a vntr
				if (is_vntr) {
					//stria_VNTR *vntr = (stria_VNTR *)PyObject_CallObject((PyObject *)&stria_VNTRType, NULL);
					stria_ETR *vntr = PyObject_New(stria_ETR, &stria_ETRType);
					vntr->motif = (char *)malloc(j + 1);
					memcpy(vntr->motif, self->seq+current_start, j);
					vntr->motif[j] = '\0';
					vntr->mlen = j;
					vntr->seqid = self->seqname;
					Py_INCREF(vntr->seqid);
					vntr->repeats = repeats;
					vntr->length = repeats * j;
					vntr->start = current_start + 1;
					vntr->end = current_start + vntr->length;
					self->next_start = vntr->end;
					return (PyObject *)vntr;
				}
			}

			i = current_start;
		}
	}

	return NULL;
}

static PyObject* stria_vntrminer_as_list(stria_VNTRMiner *self) {
	PyObject *vntrs = PyList_New(0);
	PyObject *tmp;
	Py_ssize_t current_start;
	int replen;
	int repeats;
	int length;
	int is_vntr;

	char* motif = (char *)malloc(self->max_motif + 1);

	for (Py_ssize_t i = 0; i < self->size; ++i) {
		if (self->seq[i] == 78) {
			continue;
		}

		current_start = i;
		for (int j = self->min_motif; j <= self->max_motif; ++j) {
			while (self->seq[i] == self->seq[i+j]) {
				++i;
			}

			replen = i + j - current_start;
			repeats = replen/j;

			if (repeats >= self->min_repeat) {
				//check motif is real motif with length >= min motif size
				const char *p = self->seq+current_start;
				is_vntr = 1;
				for (int k = 1; k < self->min_motif; ++k) {
					int l = 0;
					while ((p[l] == p[l+k]) && (l+k < j)) {
						++l;
					}
					if (l + k == j) {
						is_vntr = 0;
						break;
					}
				}

				if (is_vntr) {
					memcpy(motif, self->seq+current_start, j);
					motif[j] = '\0';
					length = repeats * j;
					tmp = Py_BuildValue("Onnsiii", self->seqname, current_start+1, current_start+length, motif, j, repeats, length);
					PyList_Append(vntrs, tmp);
					Py_DECREF(tmp);

					i -= replen % j;

					break;
				}
			}

			i = current_start;
		}
	}
	free(motif);
	return vntrs;
}

static PyMethodDef stria_vntrminer_methods[] = {
	{"as_list", (PyCFunction)stria_vntrminer_as_list, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

PyTypeObject stria_VNTRMinerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "VNTRMiner",                        /* tp_name */
    sizeof(stria_VNTRMiner),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)stria_vntrminer_dealloc,   /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc)stria_vntrminer_repr,                              /* tp_repr */
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
    "find minisatellites from DNA sequence",                              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    (getiterfunc)stria_vntrminer_iter,     /* tp_iter */
    (iternextfunc)stria_vntrminer_next,    /* tp_iternext */
    stria_vntrminer_methods,          /* tp_methods */
    0,          /* tp_members */
    0,                               /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    0,            /* tp_alloc */
    stria_vntrminer_new,              /* tp_new */
};
