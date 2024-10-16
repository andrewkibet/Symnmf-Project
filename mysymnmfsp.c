#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

// Function to normalize a matrix
static PyObject* norm(PyObject* self, PyObject* args) {
    PyArrayObject *X;
    int N, dim;

    if (!PyArg_ParseTuple(args, "Oii", &X, &N, &dim))
        return NULL;

    // Assuming X is a 2D array
    npy_intp* shape = PyArray_SHAPE(X);
    PyObject *norm_X = PyArray_SimpleNew(2, shape, NPY_DOUBLE);
    
    // Normalize each row
    for (int i = 0; i < N; i++) {
        double norm_val = 0;
        for (int j = 0; j < dim; j++) {
            double val = *(double*)PyArray_GETPTR2(X, i, j);
            norm_val += val * val;
        }
        norm_val = sqrt(norm_val);
        for (int j = 0; j < dim; j++) {
            double* elem = (double*)PyArray_GETPTR2(X, i, j);
            *elem /= norm_val;
        }
    }

    return norm_X;
}

// Placeholder SymNMF function
static PyObject* symnmf(PyObject* self, PyObject* args) {
    int N, dim, k;

    if (!PyArg_ParseTuple(args, "iii", &N, &dim, &k))
        return NULL;

    // Placeholder: returns a random matrix for testing purposes
    npy_intp dims[2] = {N, k};
    PyObject *result = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            double rand_val = rand() / (double)RAND_MAX;
            *(double*)PyArray_GETPTR2(result, i, j) = rand_val;
        }
    }

    return result;
}

// Method definitions
static PyMethodDef SymNMFMethods[] = {
    {"norm", norm, METH_VARARGS, "Normalize the matrix."},
    {"symnmf", symnmf, METH_VARARGS, "Symmetric Non-negative Matrix Factorization."},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "mysymnmfsp",
    NULL,
    -1,
    SymNMFMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_mysymnmfsp(void) {
    import_array();  // Initialize NumPy
    return PyModule_Create(&symnmfmodule);
}
