#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "etr.h"
#include "ssr.h"
#include "itr.h"
#include "math.h"
#include "vntr.h"
#include "version.h"

PyObject *test(PyObject *self, PyObject *args, PyObject *kwargs) {
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
	Py_ssize_t bs[7];

	for (j = 0; j < 7; ++j) {
		bs[j] = size - j;
	}

	for (i = 0; i < size; ++i) {
		if (seq[i] == 78) {
			continue;
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

			while ((i < bs[j]) && (seq[i] == seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= min_lens[j]) {
				memcpy(motif, seq+current_start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;
				ssr_end = current_start+length;
				tmp = Py_BuildValue("snnsiii", name, current_start+1, ssr_end, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = ssr_end;
				break;
			} else {
				i = current_start;
			}
		}
	}

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
	int h;
	int replen;
	int repeats;
	int length;
	char motif[7];
	int min_lens[] = {0, 12, 14, 15, 16, 20, 24};
	int hash[10] = {0};

	static char* keywords[] = {"name", "seq", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", keywords, &name, &seqobj)) {
		return NULL;
	}

	Py_INCREF(seqobj);
	seq = PyUnicode_AsUTF8AndSize(seqobj, &size);

	//boundaries for each motif length
	Py_ssize_t bs[7];

	for (j = 0; j < 7; ++j) {
		bs[j] = size - j;
	}

	for (i = 0; i < size; ++i) {
		if (seq[i] == 78) {
			continue;
		}

		/*current_start = i;
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
		}*/

		for (j=1; j < 7; ++j) {

			//find candidate motif
			hash[j] = hash[j]
			


			while ((i < bs[j]) && (seq[i] == seq[i+j])) {
				++i;
			}

			replen = i + j - current_start;

			if (replen >= min_lens[j]) {
				memcpy(motif, seq+current_start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;
				ssr_end = current_start+length;
				tmp = Py_BuildValue("snnsiii", name, current_start+1, ssr_end, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i = ssr_end;
				break;
			} else {
				i = current_start;
			}
		}
	}

	Py_DECREF(seqobj);
	return ssrs;
}

PyObject *version(PyObject *self) {
	return PyUnicode_FromString(STRIA_VERSION);
}

static PyMethodDef module_methods[] = {
	{"test", (PyCFunction)test, METH_VARARGS | METH_KEYWORDS, NULL},
	{"version", (PyCFunction)version, METH_NOARGS, NULL},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_stria = {
	PyModuleDef_HEAD_INIT,
	"stria",
	"A python package for short tandem repeat identification and analysis",
	-1,
	module_methods,
};

static PyObject *stria_module_init(void) {
	PyObject *module;
	module = PyModule_Create(&module_stria);

	if (module == NULL) {
		return NULL;
	}

	//TRE
	if (PyType_Ready(&stria_ETRType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stria_ETRType);
	PyModule_AddObject(module, "ETR", (PyObject *)&stria_ETRType);

	//SSR
	if (PyType_Ready(&stria_SSRMinerType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stria_SSRMinerType);
	PyModule_AddObject(module, "SSRMiner", (PyObject *)&stria_SSRMinerType);

	//VNTR
	if (PyType_Ready(&stria_VNTRMinerType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stria_VNTRMinerType);
	PyModule_AddObject(module, "VNTRMiner", (PyObject *)&stria_VNTRMinerType);

	//ITR
	if (PyType_Ready(&stria_ITRType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stria_ITRType);
	PyModule_AddObject(module, "ITR", (PyObject *)&stria_ITRType);

	if (PyType_Ready(&stria_ITRMinerType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stria_ITRMinerType);
	PyModule_AddObject(module, "ITRMiner", (PyObject *)&stria_ITRMinerType);

	return module;
}

PyMODINIT_FUNC PyInit_stria() {
	return stria_module_init();
}
