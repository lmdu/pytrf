#include <Python.h>
#include "ssr.h"

PyObject *test(PyObject *self, PyObject *args, PyObject *kwargs) {
	return Py_BuildValue("s", "hello");
}

static PyMethodDef module_methods[] = {
	{"test", (PyCFunction)test, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_stripy = {
	PyModuleDef_HEAD_INIT,
	"stripy",
	"A python package for short tandem repeat identification",
	-1,
	module_methods,
};

static PyObject *strit_module_init(void) {
	PyObject *module;
	module = PyModule_Create(&module_stripy);

	if (module == NULL) {
		return NULL;
	}

	if (PyType_Ready(&stripy_SSRMinerType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stripy_SSRMinerType);
	PyModule_AddObject(module, "SSRMiner", (PyObject *)&stripy_SSRMinerType);

	if (PyType_Ready(&stripy_SSRType) < 0) {
		return NULL;
	}
	Py_INCREF((PyObject *)&stripy_SSRType);
	PyModule_AddObject(module, "SSR", (PyObject *)&stripy_SSRType);

	return module;
}

PyMODINIT_FUNC PyInit_stripy() {
	return strit_module_init();
}
