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

static PyObject* npy_i64(const std::int64_t* data, npy_intp rows, npy_intp cols)
{
    npy_intp dims[2] = { rows, cols };
    PyObject* arr = PyArray_SimpleNew(2, dims, NPY_INT64);
    if (!arr) return NULL;
    std::memcpy(PyArray_DATA((PyArrayObject*)arr), data, (size_t)(rows * cols) * sizeof(std::int64_t));
    return arr;
}

static PyObject* npy_f64_1d(const double* data, npy_intp n)
{
    PyObject* arr = PyArray_SimpleNew(1, &n, NPY_FLOAT64);
    if (!arr) return NULL;
    std::memcpy(PyArray_DATA((PyArrayObject*)arr), data, (size_t)n * sizeof(double));
    return arr;
}

static PyObject* npy_i64_1d(const std::int64_t* data, npy_intp n)
{
    PyObject* arr = PyArray_SimpleNew(1, &n, NPY_INT64);
    if (!arr) return NULL;
    std::memcpy(PyArray_DATA((PyArrayObject*)arr), data, (size_t)n * sizeof(std::int64_t));
    return arr;
}

static PyObject* result_to_dict(TDResult& R)
{
    if (!R.error.empty()) {
        PyErr_SetString(PyExc_RuntimeError, R.error.c_str());
        return NULL;
    }

    PyObject* splits = PyList_New((Py_ssize_t)R.triangle_splits.size());
    if (!splits) return NULL;
    for (size_t i = 0; i < R.triangle_splits.size(); ++i) {
        const TDSplit& s = R.triangle_splits[i];
        PyObject* d = Py_BuildValue(
            "{s:i,s:i,s:i,s:i,s:N,s:N,s:N,s:N}",
            "split_id", s.split_id,
            "iteration", s.iteration,
            "priority", s.priority,
            "is_lookahead", s.is_lookahead,
            "parent", npy_f64(&s.parent.v[0][0], 3, 2),
            "child1", npy_f64(&s.child1.v[0][0], 3, 2),
            "child2", npy_f64(&s.child2.v[0][0], 3, 2),
            "split_line", npy_f64(&s.split_line.a[0], 2, 2));
        if (!d) { Py_DECREF(splits); return NULL; }
        PyList_SET_ITEM(splits, (Py_ssize_t)i, d);
    }

    PyObject* lines = PyList_New((Py_ssize_t)R.line_splits.size());
    if (!lines) { Py_DECREF(splits); return NULL; }
    for (size_t i = 0; i < R.line_splits.size(); ++i) {
        const TDLineSplit& s = R.line_splits[i];
        std::int64_t qp[2] = { s.qx, s.qy };
        PyObject* d = Py_BuildValue(
            "{s:i,s:i,s:i,s:N,s:N,s:N,s:N,s:N}",
            "split_id", s.split_id,
            "triangle_split_id", s.triangle_split_id,
            "iteration", s.iteration,
            "original_edge", npy_f64(&s.original_edge.a[0], 2, 2),
            "segment1", npy_f64(&s.segment1.a[0], 2, 2),
            "segment2", npy_f64(&s.segment2.a[0], 2, 2),
            "split_point", npy_f64_1d(s.split_point, 2),
            "quantized_point", npy_i64_1d(qp, 2));
        if (!d) { Py_DECREF(splits); Py_DECREF(lines); return NULL; }
        PyList_SET_ITEM(lines, (Py_ssize_t)i, d);
    }

    PyObject* marks = PyList_New((Py_ssize_t)R.marks.size());
    if (!marks) { Py_DECREF(splits); Py_DECREF(lines); return NULL; }
    for (size_t i = 0; i < R.marks.size(); ++i) {
        const TDMark& m = R.marks[i];
        PyObject* sides = PyList_New(0);
        PyObject* d = Py_BuildValue(
            "{s:i,s:i,s:s,s:N,s:d,s:N,s:N}",
            "leaf_id", m.leaf_id,
            "region_code", m.region_code,
            "region", m.region.c_str(),
            "vertices", npy_f64(&m.vertices.v[0][0], 3, 2),
            "area", m.area,
            "quantized", npy_i64(&m.qv[0][0], 3, 2),
            "polygon_sides", sides);
        if (!d) { Py_DECREF(splits); Py_DECREF(lines); Py_DECREF(marks); return NULL; }
        PyList_SET_ITEM(marks, (Py_ssize_t)i, d);
    }

    PyObject* walk = PyList_New((Py_ssize_t)R.walk.size());
    if (!walk) { Py_DECREF(splits); Py_DECREF(lines); Py_DECREF(marks); return NULL; }
    for (size_t i = 0; i < R.walk.size(); ++i) {
        const TDWalkSeg& w = R.walk[i];
        double exact[4] = { w.sx, w.sy, w.tx, w.ty };
        std::int64_t q[4] = { w.qsx, w.qsy, w.qtx, w.qty };
        PyObject* d = Py_BuildValue(
            "{s:i,s:N,s:N}",
            "seq", w.seq,
            "exact", npy_f64(exact, 2, 2),
            "quantized", npy_i64(q, 2, 2));
        if (!d) {
            Py_DECREF(splits); Py_DECREF(lines); Py_DECREF(marks); Py_DECREF(walk);
            return NULL;
        }
        PyList_SET_ITEM(walk, (Py_ssize_t)i, d);
    }

    PyObject* stats = Py_BuildValue(
        "{s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:d,s:d,s:O,s:O,s:s}",
        "n_polygon", R.stats.n_polygon,
        "n_splits", R.stats.n_splits,
        "n_line_splits", R.stats.n_line_splits,
        "n_leaves", R.stats.n_leaves,
        "n_interior", R.stats.n_interior,
        "n_exterior", R.stats.n_exterior,
        "perwalk_seeds", R.stats.perwalk_seeds,
        "flood_added", R.stats.flood_added,
        "flood_unreached", R.stats.flood_unreached,
        "inside_area", R.stats.inside_area,
        "polygon_area", R.stats.polygon_area,
        "area_ok", R.stats.area_ok ? Py_True : Py_False,
        "walkable", R.stats.walkable ? Py_True : Py_False,
        "engine", "cpp");
    if (!stats) {
        Py_DECREF(splits); Py_DECREF(lines); Py_DECREF(marks); Py_DECREF(walk);
        return NULL;
    }

    PyObject* out = Py_BuildValue(
        "{s:N,s:N,s:N,s:N,s:N}",
        "triangle_splits", splits,
        "line_splits", lines,
        "marks", marks,
        "walk", walk,
        "stats", stats);
    return out;
}

static PyObject* py_tridecomp(PyObject* /*self*/, PyObject* args, PyObject* kwargs)
{
    static const char* kwlist[] = { "p", "t", "verbose", NULL };
    PyObject *pobj = NULL, *tobj = NULL;
    int verbose = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|p:tridecomp",
                                     (char**)kwlist, &pobj, &tobj, &verbose))
        return NULL;

    PyArrayObject* pa = (PyArrayObject*)PyArray_FROM_OTF(pobj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject* ta = (PyArrayObject*)PyArray_FROM_OTF(tobj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (!pa || !ta) {
        Py_XDECREF(pa); Py_XDECREF(ta);
        PyErr_SetString(PyExc_TypeError, "p and t must be array-like float64");
        return NULL;
    }
    if (PyArray_NDIM(pa) != 2 || PyArray_DIM(pa, 1) != 2) {
        Py_DECREF(pa); Py_DECREF(ta);
        PyErr_SetString(PyExc_ValueError, "p must have shape (N, 2)");
        return NULL;
    }
    if (PyArray_NDIM(ta) != 2 || PyArray_DIM(ta, 0) != 3 || PyArray_DIM(ta, 1) != 2) {
        Py_DECREF(pa); Py_DECREF(ta);
        PyErr_SetString(PyExc_ValueError, "t must have shape (3, 2)");
        return NULL;
    }

    const int n = (int)PyArray_DIM(pa, 0);
    const double* pxy = (const double*)PyArray_DATA(pa);
    const double* txy = (const double*)PyArray_DATA(ta);

    // Copies if the array is not C-contiguous row-major — FROM_OTF IN_ARRAY
    // already made a C-contiguous copy if needed.
    TDResult R = tridecomp_from_xy(pxy, n, txy, verbose != 0);
    Py_DECREF(pa);
    Py_DECREF(ta);
    return result_to_dict(R);
}

static PyMethodDef methods[] = {
    { "tridecomp", (PyCFunction)py_tridecomp, METH_VARARGS | METH_KEYWORDS,
      "tridecomp(p, t, verbose=False) — EPEC C++ engine" },
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "tridecomp._cpp",
    "CGAL EPEC triangular decomposition",
    -1,
    methods
};

PyMODINIT_FUNC PyInit__cpp(void)
{
    import_array();
    if (PyErr_Occurred()) return NULL;
    return PyModule_Create(&moduledef);
}
