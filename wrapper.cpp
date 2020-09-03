
#include <Python.h>
#include <stdio.h>
#include "sharedmatting.h"

// https://docs.scipy.org/doc/numpy/reference/c-api/types-and-structures.html
typedef struct {
    PyObject_HEAD
    uint8_t *data;
    int nd;
    int64_t *dimensions;
    int64_t *strides;
} PyArrayObject;

// https://python3-cookbook.readthedocs.io/zh_CN/latest/chapters/p15_c_extensions.html
static PyObject* py_solve_alpha(PyObject *self, PyObject *args)
{
    PyArrayObject *im, *trimap, *result;
    if (!PyArg_ParseTuple(args, "OOO", &im, &trimap, &result))
        return NULL;

    SharedMatting sm;
    sm.loadImage(im->data, trimap->data, im->dimensions[1], im->dimensions[0]);
    sm.solveAlpha(result->data);

    Py_RETURN_NONE;
}

/* Module method table */
static PyMethodDef SampleMethods[] = {
  {"solve_alpha", py_solve_alpha, METH_VARARGS, "solve_alpha"},
  { NULL, NULL, 0, NULL}
};

/* Module structure */
static struct PyModuleDef samplemodule = {
  PyModuleDef_HEAD_INIT,
  "sample",           /* name of module */
  "A sample module",  /* Doc string (may be NULL) */
  -1,                 /* Size of per-interpreter state or -1 */
  SampleMethods       /* Method table */
};

/* Module initialization function */
PyMODINIT_FUNC
PyInit_sharedmatting(void) {
    return PyModule_Create(&samplemodule);
}
