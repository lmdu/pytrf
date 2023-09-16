#ifndef PYTRF_COMPAT_H
#define PYTRF_COMPAT_H

#include <Python.h>

#ifndef _Py_CAST
#define _Py_CAST(type, expr) ((type)(expr))
#endif

#ifndef _PyObject_CAST
#define _PyObject_CAST(op) _Py_CAST(PyObject*, op)
#endif

#if PY_VERSION_HEX < 0x030A00A3 && !defined(Py_NewRef)
static inline PyObject* _Py_NewRef(PyObject *obj) {
	Py_INCREF(obj);
	return obj;
}
#define Py_NewRef(obj) _Py_NewRef(_PyObject_CAST(obj))
#endif

#if PY_VERSION_HEX < 0x030A00A3 && !defined(Py_XNewRef)
static inline PyObject* _Py_XNewRef(PyObject *obj) {
	Py_XINCREF(obj);
	return obj;
}
#define Py_XNewRef(obj) _Py_XNewRef(_PyObject_CAST(obj))
#endif

#endif
