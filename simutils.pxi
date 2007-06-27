cimport c_python
cimport c_numpy
# Numpy must be initialized
c_numpy.import_array()

import sys
import time
from Waitbar import Waitbar
import sim

cdef double* DoubleData(c_numpy.ndarray M):
    return <double *>M.data

cdef char* CharData(c_numpy.ndarray M):
    return <char *>M.data
    
cdef int* IntData(c_numpy.ndarray M):
    return <int *>M.data
    

cdef int Dim0(c_numpy.ndarray M):
    return M.dimensions[0]

cdef int Dim1(c_numpy.ndarray M):
    return M.dimensions[1]

cdef extern from "randmtzig.c":
    void init_by_int(int)
    void init_by_entropy()

    double randu()
    double randn()
    double rande()

cdef extern from "math.h":
    double floor(double)
    double exp(double)
    double log(double)
    double tanh(double)


