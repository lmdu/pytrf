#include <Python.h>
#include "tandem.h"

static PyMethodDef module_methods[] = {
	{"find_ssrs", (PyCFunction)find_short_tandem_repeats, METH_VARARGS | METH_KEYWORDS, NULL},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_strit = {
	PyModuleDef_HEAD_INIT,
	"strit",
	"A python C extension for finding tandem repeats",
	-1,
	module_methods,
};

static PyObject *strit_module_init(void) {
	PyObject *module;
	module = PyModule_Create(&module_strit);

	if (module == NULL) {
		return NULL;
	}

	return module;
}

PyMODINIT_FUNC PyInit_strit() {
	return strit_module_init();
}
