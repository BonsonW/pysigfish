#cython: language_level=3
from libc.stdio cimport *
from libc.stdint cimport *
from libc.stdlib cimport *

cdef extern from "sigfish.h":


    # init sigfish
    int initsigfish();
    