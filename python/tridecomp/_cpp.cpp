#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "tridecomp_api.h"

#include <cstring>
#include <vector>

static PyObject* npy_f64(const double* data, npy_intp rows, npy_intp cols)
{
    npy_intp dims[2] = { rows, cols };
    PyObject* arr = PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (!arr) return NULL;
    std::memcpy(PyArray_DATA((PyArrayObject*)arr), data, (size_t)(rows * cols) * sizeof(double));
    return arr;
}

static PyObject* npy_f64_1d(const double* data, npy_intp n)
{
    PyObject* arr = PyArray_SimpleNew(1, &n, NPY_FLOAT64);
    if (!arr) return NULL;
    std::memcpy(PyArray_DATA((PyArrayObject*)arr), data, (size_t)n * sizeof(double));
    return arr;
}

static PyObject* line_split_dict(const TDLineSplit& s)
{
    return Py_BuildValue(
        "{s:N,s:N,s:(NN)}",
        "edge", npy_f64(&s.original_edge.a[0], 2, 2),
        "point", npy_f64_1d(s.split_point, 2),
        "seg", npy_f64(&s.segment1.a[0], 2, 2), npy_f64(&s.segment2.a[0], 2, 2));
}

static PyObject* events_to_list(TDTessellation& R)
{
    if (!R.error.empty()) {
        PyErr_SetString(PyExc_RuntimeError, R.error.c_str());
        return NULL;
    }

    PyObject* events = PyList_New((Py_ssize_t)R.events.size());
    if (!events) return NULL;
    for (size_t i = 0; i < R.events.size(); ++i) {
        const TDSplitEvent& e = R.events[i];
        PyObject* lines = PyList_New((Py_ssize_t)e.line_splits.size());
        if (!lines) { Py_DECREF(events); return NULL; }
        for (size_t k = 0; k < e.line_splits.size(); ++k) {
            PyObject* ls = line_split_dict(e.line_splits[k]);
            if (!ls) { Py_DECREF(lines); Py_DECREF(events); return NULL; }
            PyList_SET_ITEM(lines, (Py_ssize_t)k, ls);
        }
        PyObject* d = Py_BuildValue(
            "{s:s,s:s,s:s,s:N,s:N,s:(NN),s:N}",
            "id", e.id.c_str(),
            "child0", e.child0.c_str(),
            "child1", e.child1.c_str(),
            "parent", npy_f64(&e.parent.v[0][0], 3, 2),
            "cut", npy_f64(&e.cut.a[0], 2, 2),
            "child", npy_f64(&e.child[0].v[0][0], 3, 2), npy_f64(&e.child[1].v[0][0], 3, 2),
            "line_splits", lines);
        if (!d) { Py_DECREF(events); return NULL; }
        PyList_SET_ITEM(events, (Py_ssize_t)i, d);
    }
    return events;
}

static int as_segments(PyObject* obj, std::vector<double>& flat, int* nseg)
{
    PyArrayObject* a = (PyArrayObject*)PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (!a) return -1;
    const int nd = PyArray_NDIM(a);
    const npy_intp* dims = PyArray_DIMS(a);
    if (nd == 3 && dims[1] == 2 && dims[2] == 2) {
        *nseg = (int)dims[0];
    } else if (nd == 2 && dims[1] == 4) {
        *nseg = (int)dims[0];
    } else {
        Py_DECREF(a);
        PyErr_SetString(PyExc_ValueError, "segments must have shape (M, 2, 2) or (M, 4)");
        return -1;
    }
    const npy_intp n = (npy_intp)(*nseg) * 4;
    flat.resize((size_t)n);
    std::memcpy(flat.data(), PyArray_DATA(a), (size_t)n * sizeof(double));
    Py_DECREF(a);
    return 0;
}

static PyObject* py_tessellate(PyObject* /*self*/, PyObject* args, PyObject* kwargs)
{
    static const char* kwlist[] = { "segments", "t", "start", "verbose", NULL };
    PyObject *sobj = NULL, *tobj = NULL;
    const char* start = "";
    int verbose = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|$sp:tessellate",
                                     (char**)kwlist, &sobj, &tobj, &start, &verbose))
        return NULL;

    std::vector<double> segflat;
    int nseg = 0;
    if (as_segments(sobj, segflat, &nseg) != 0)
        return NULL;

    PyArrayObject* ta = (PyArrayObject*)PyArray_FROM_OTF(tobj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (!ta) {
        PyErr_SetString(PyExc_TypeError, "t must be array-like float64");
        return NULL;
    }
    if (PyArray_NDIM(ta) != 2 || PyArray_DIM(ta, 0) != 3 || PyArray_DIM(ta, 1) != 2) {
        Py_DECREF(ta);
        PyErr_SetString(PyExc_ValueError, "t must have shape (3, 2)");
        return NULL;
    }
    const double* txy = (const double*)PyArray_DATA(ta);
    const double* sxy = nseg ? segflat.data() : nullptr;
    TDTessellation R = tessellate_from_xy(sxy, nseg, txy, start ? start : "", verbose != 0);
    Py_DECREF(ta);
    return events_to_list(R);
}

static PyMethodDef methods[] = {
    { "tessellate", (PyCFunction)py_tessellate, METH_VARARGS | METH_KEYWORDS,
      "tessellate(segments, t, start='', verbose=False) — EPEC C++ engine" },
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "tridecomp._cpp",
    "CGAL EPEC tessellate(segments, t)",
    -1,
    methods
};

PyMODINIT_FUNC PyInit__cpp(void)
{
    import_array();
    if (PyErr_Occurred()) return NULL;
    return PyModule_Create(&moduledef);
}
