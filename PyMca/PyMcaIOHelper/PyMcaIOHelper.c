#include "Python.h"
#include <./numpy/arrayobject.h>
#include <stdio.h>

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

/* Function declarations */

static PyObject *PyMcaIOHelper_fillSupaVisio(PyObject *dummy, PyObject *args);
static PyObject *PyMcaIOHelper_readAifira(PyObject *dummy, PyObject *args);

/* Functions */

static PyObject *
PyMcaIOHelper_fillSupaVisio(PyObject *self, PyObject *args)
{
    PyObject *input;
    PyArrayObject *inputArray, *outputArray;
    int nChannels = 2048;
    unsigned short *dataPointer, x, y, ch;
    int i;
    npy_intp dimensions[3];
    int maxy, maxch;
    unsigned int *outputPointer;

    if (!PyArg_ParseTuple(args, "O|d", &input, &nChannels))
        return NULL;
    inputArray = (PyArrayObject *)
                PyArray_ContiguousFromObject(input, PyArray_USHORT, 2,2);
    if (inputArray == NULL)
    {
	    struct module_state *st = GETSTATE(self);
		PyErr_SetString(st->error, "Cannot parse input array");
        return NULL;
    }
    
    dataPointer = (unsigned short *) inputArray->data;
    dataPointer++;
    dimensions[1] = *dataPointer++;
    dimensions[0] = *dataPointer++;
    dimensions[2] = nChannels;
    outputArray = (PyArrayObject *) PyArray_SimpleNew(3, dimensions, PyArray_UINT);
    PyArray_FILLWBYTE(outputArray, 0);
    /* Do the job */
    maxy=maxch=0;
    outputPointer = (unsigned int *) outputArray->data; 
    for (i=3; i<inputArray->dimensions[0]; i++)
    {
        y = *dataPointer++;
    x = *dataPointer++;
    ch = *dataPointer++;
    if (ch >= nChannels)
    {
        printf("bad reading %d\n", ch);
        continue;
    }
    *(outputPointer+ (dimensions[1] * x + y) * nChannels + ch) += 1;
    }
    Py_DECREF(inputArray);
    return PyArray_Return(outputArray);
}

static PyObject *
PyMcaIOHelper_readAifira(PyObject *self, PyObject *args)
{
    PyObject *inputFileDescriptor;
    FILE *fd;
#if PY_MAJOR_VERSION >= 3
	int fh;
#endif
    PyArrayObject *outputArray;
    int nChannels = 2048;
    unsigned short channel;
    unsigned char x, y; 
    npy_intp dimensions[3];
    unsigned int *outputPointer;
	struct module_state *st = GETSTATE(self);

    if (!PyArg_ParseTuple(args, "O", &inputFileDescriptor))
    {
		PyErr_SetString(st->error, "Error parsing input arguments");
        return NULL;
    }
#if PY_MAJOR_VERSION >= 3
    fh = PyObject_AsFileDescriptor(inputFileDescriptor);
	if (fh < 0)
    {
        return NULL;
    }
	fd = fdopen(fh, "r");
#else
	if (!PyFile_Check(inputFileDescriptor))
    {
        PyErr_SetString(st->error, "Input is not a python file descriptor object");
        return NULL;
    }
    fd = PyFile_AsFile(inputFileDescriptor);
#endif
	if (!fd)
	{
        PyErr_SetString(st->error, "Cannot obtain FILE* from object");
        return NULL;
	}
    dimensions[0] = 128;
    dimensions[1] = 128;
    dimensions[2] = nChannels;

    outputArray = (PyArrayObject *) PyArray_SimpleNew(3, dimensions, PyArray_UINT);
    PyArray_FILLWBYTE(outputArray, 0);

    /* Do the job */
    outputPointer = (unsigned int *) outputArray->data;
    while(fscanf(fd, "%2c%c%c", (char *)&channel, &x, &y) == 3)
    {
        if (channel >= nChannels)
        {
            printf("bad reading %d\n", channel);
            continue;
        }
        if (x < 1)
            continue;
        if (y < 1)
            continue;
        if (x > 128)
        {
            printf("bad X reading %d\n", x);
            break;
            continue;
        }
        if (y > 128)
        {
            printf("bad Y reading %d\n", y);
            break;
            continue;
        }
        x -= 1;
        y -= 1;
        /* normally pixe data are in the second channel */
        if (channel > 1023)
        {
            channel -= 1024;
        }
        else
            channel += 1024;

        *(outputPointer + (dimensions[1] * x + y) * nChannels + channel) += 1;
    }
    return PyArray_Return(outputArray);
}

/* Module methods */

static PyMethodDef PyMcaIOHelper_methods[] = {
    {"fillSupaVisio", PyMcaIOHelper_fillSupaVisio, METH_VARARGS},
    {"readAifira", PyMcaIOHelper_readAifira, METH_VARARGS},
	{NULL, NULL}
};

/* ------------------------------------------------------- */


/* Module initialization */

#if PY_MAJOR_VERSION >= 3

static int PyMcaIOHelper_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int PyMcaIOHelper_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "PyMcaIOHelper",
        NULL,
        sizeof(struct module_state),
        PyMcaIOHelper_methods,
        NULL,
        PyMcaIOHelper_traverse,
        PyMcaIOHelper_clear,
        NULL
};

#define INITERROR return NULL

PyObject *
PyInit_PyMcaIOHelper(void)

#else
#define INITERROR return

void
initPyMcaIOHelper(void)
#endif
{
	struct module_state *st;
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("PyMcaIOHelper", PyMcaIOHelper_methods);
#endif

    if (module == NULL)
        INITERROR;
    st = GETSTATE(module);

    st->error = PyErr_NewException("PyMcaIOHelper.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    import_array();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
