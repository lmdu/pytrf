#include <Python.h>
#include "etr.h"
#include "ssr.h"
#include "itr.h"
#include "vntr.h"
#include "version.h"

PyObject *test(PyObject *self, PyObject *args, PyObject *kwargs) {
	PyObject *ssrs = PyList_New(0);
	PyObject *tmp;
	PyObject *seqobj;
	PyObject *seqname;

	Py_ssize_t start;
	Py_ssize_t size;

	int replen;
	int repeats;
	int length;
	char motif[7];
	int min_lens[] = {0, 12, 14, 15, 16, 20, 24};

	static char* keywords[] = {"name", "seq", NULL};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", keywords, &seqname, &seqobj)) {
		return NULL;
	}

	size = PyObject_Length(seqobj);

	Py_UCS1* seq = PyUnicode_1BYTE_DATA(seqobj);

	for (Py_ssize_t i = 0; i < size; ++i) {
		if (seq[i] == 78) {
			continue;
		}

		start = i;
		for (int j = 1; j < 7; ++j) {
			while (seq[i] == seq[i+j]) {
				++i;
			}

			replen = i + j - start;

			if (replen >= min_lens[j]) {
				memcpy(motif, seq+start, j);
				motif[j] = '\0';
				repeats = replen/j;
				length = repeats * j;

				tmp = Py_BuildValue("Onnsiii", seqname, start+1, start+length, motif, j, repeats, length);
				PyList_Append(ssrs, tmp);
				Py_DECREF(tmp);

				i -= replen % j;

				break;
			} else {
				i = start;
			}
		}
	}

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
