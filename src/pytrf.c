#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "atr.h"
#include "etr.h"
#include "str.h"
#include "itr.h"
#include "gtr.h"
#include "math.h"
#include "version.h"

int nt_code[128] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

PyObject *test(PyObject *self, PyObject *args, PyObject *kwargs) {
	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;

	PyObject *name;
	PyObject *seqobj;

	const char *seq;

	//current start position
	Py_ssize_t cur;

	//ssr end
	Py_ssize_t se;

	Py_ssize_t size;
	Py_ssize_t i;

	int j;

	//repeat length
	int rl;

	//repeat number
	int r;

	//ssr length
	int sl;

	//motif
	char motif[7];

	//min length for each repeat type
	int ml[] = {0, 12, 14, 15, 16, 20, 24};

	//boundaries for each motif length
	Py_ssize_t bs[7];

	static char* keywords[] = {"name", "seq", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", keywords, &name, &seqobj)) {
		return NULL;
	}

	Py_INCREF(name);
	Py_INCREF(seqobj);
	seq = PyUnicode_AsUTF8AndSize(seqobj, &size);

	for (j = 0; j < 7; ++j) {
		bs[j] = size - j;
	}

	for (i = 0; i < size; ++i) {
		//pass unkown bases
		if (seq[i] == 'N' || seq[i] == 'n') {
			continue;
		}

		cur = i;

		j = 1;

		/*while ((i < bs[1]) && (seq[i] == seq[i+1])) {
			++i;
		}

		rl = i + 1 - cur;

		if (rl < ml[j]) {
			j = rl + 1;
			i = cur;

			if (j > 6) {
				i += rl - 6;
				continue;
			}
		}*/

		for (; j < 7; ++j) {
			while ((i < bs[j]) && (seq[i] == seq[i+j])) {
				++i;
			}

			if (i == cur) {
				continue;
			}

			rl = i + j - cur;

			if (rl >= ml[j]) {
				memcpy(motif, seq+cur, j);
				motif[j] = '\0';
				r = rl / j;
				sl = r * j;
				se = cur + sl;
				tmp = Py_BuildValue("Onnsiii", name, cur+1, se, motif, j, r, sl);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = se;
				break;
			} else {
				i = cur;
			}
		}
	}

	Py_DECREF(name);
	Py_DECREF(seqobj);
	return ssrs;
}

PyObject *test_hash(PyObject *self, PyObject *args, PyObject *kwargs) {
	PyObject *ssrs = PyList_New(0);
	PyObject *name;
	PyObject *tmp;
	PyObject *seqobj;

	const char *seq;
	Py_ssize_t current_start;
	Py_ssize_t ssr_end;
	Py_ssize_t size;
	Py_ssize_t i;

	int j;
	int k;
	int replen;
	int repeats;
	int length;
	char motif[7];
	int min_lens[] = {0, 12, 14, 15, 16, 20, 24};

	static char* keywords[] = {"name", "seq", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", keywords, &name, &seqobj)) {
		return NULL;
	}

	Py_INCREF(seqobj);
	seq = PyUnicode_AsUTF8AndSize(seqobj, &size);

	//boundaries for each motif length
	Py_ssize_t bs[7] = {0};

	for (j = 1; j < 7; ++j) {
		bs[j] = size - j;
	}

	Py_ssize_t hash = 0;

	int ps[7] = {0, 0, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6};

	//split
	int ss[7] = {0, 0, 1e+8, 1e+6, 1e+4, 1e+2, 1};

	//max pow
	Py_ssize_t mp = 1e+11;

	for (i = 0; i < size; ++i) {
		if (seq[i] == 78) {
			continue;
		}

		//rolling hash
		if (i > bs[6]) {
			break;
		}
		
		if (!hash) {
			for (k = 0; k < 12; ++k) {
				hash = hash*10 + nt_code[(unsigned char)seq[i+k]];
			}
		} else {
			hash = (hash % mp) * 10 + nt_code[(unsigned char)seq[i]];
		}

		current_start = i;
		j = 1;

		while ((i < bs[1]) && (seq[i] == seq[i+1])) {
			++i;
		}

		replen = i - current_start + 1;

		if (replen < min_lens[j]) {
			j = replen + 1;
			i = current_start;

			if (j > 6) {
				i += replen - 6;
				continue;
			}
		}

		for (; j < 7; ++j) {

			//find candidate motif
			/*if (j > 1) {
				Py_ssize_t num = hash/ss[j];

				if (num/ps[j] != num%ps[j]) {
					continue;
				}
			}*/

			while ((i < bs[j]) && (seq[i] == seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= min_lens[j]) {
				//printf("%zd\n", hash);
				memcpy(motif, seq+current_start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;
				ssr_end = current_start+length;
				tmp = Py_BuildValue("snnsiii", name, current_start+1, ssr_end, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = ssr_end;
				hash = 0;
				break;
			} else {
				i = current_start;
			}
		}
	}

	Py_DECREF(seqobj);
	return ssrs;
}

/*PyObject *version(PyObject *self) {
	return PyUnicode_FromString(STRIA_VERSION);
}*/

static PyMethodDef pytrf_methods[] = {
	{"test", (PyCFunction)test, METH_VARARGS | METH_KEYWORDS, NULL},
	{"test_hash", (PyCFunction)test_hash, METH_VARARGS | METH_KEYWORDS, NULL},
	//{"version", (PyCFunction)version, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef pytrf_module = {
	PyModuleDef_HEAD_INIT,
	"pytrf",
	"A python package for finding tandem repeat sequence",
	-1,
	pytrf_methods,
};

static PyObject *pytrf_module_init(void) {
	PyObject *module;
	module = PyModule_Create(&pytrf_module);

	if (module == NULL) {
		return NULL;
	}

	//version
	PyModule_AddStringConstant(module, "__version__", PYTRF_VERSION);

	//ETR
	if (PyType_Ready(&pytrf_ETRType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&pytrf_ETRType);
	PyModule_AddObject(module, "ETR", (PyObject *)&pytrf_ETRType);

	//ATR
	if (PyType_Ready(&pytrf_ATRType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&pytrf_ATRType);
	PyModule_AddObject(module, "ATR", (PyObject *)&pytrf_ATRType);

	//STR
	if (PyType_Ready(&pytrf_STRFinderType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&pytrf_STRFinderType);
	PyModule_AddObject(module, "STRFinder", (PyObject *)&pytrf_STRFinderType);

	//GTR
	if (PyType_Ready(&pytrf_GTRFinderType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&pytrf_GTRFinderType);
	PyModule_AddObject(module, "GTRFinder", (PyObject *)&pytrf_GTRFinderType);

	//ITR
	if (PyType_Ready(&pytrf_ITRFinderType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&pytrf_ITRFinderType);
	PyModule_AddObject(module, "ITRFinder", (PyObject *)&pytrf_ITRFinderType);

	return module;
}

PyMODINIT_FUNC PyInit_pytrf() {
	return pytrf_module_init();
}
