#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "atr.h"
#include "etr.h"
#include "str.h"
#include "itr.h"
#include "gtr.h"
#include "math.h"
#include "version.h"

PyObject *pytrf_version(PyObject *self) {
	return PyUnicode_FromString(PYTRF_VERSION);
}

static PyMethodDef pytrf_methods[] = {
	{"version", (PyCFunction)pytrf_version, METH_NOARGS, NULL},
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
	PyModule_AddObject(module, "ATRFinder", (PyObject *)&pytrf_ITRFinderType);

	return module;
}

PyMODINIT_FUNC PyInit_pytrf() {
	return pytrf_module_init();
}
