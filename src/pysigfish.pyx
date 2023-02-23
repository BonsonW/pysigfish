# distutils: language = c
# cython: language_level=3
# cython: profile=True
import sys
import time
import logging
import copy
from libc.stdlib cimport malloc, free
from libc.string cimport strdup
cimport sigfish
# Import the Python-level symbols of numpy
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()


cdef class start:
    '''
    Creates a new sigfish object

    '''
    def __cinit__(self, threads=8, DEBUG=0):
        '''
        C init
        '''

        # sets up logging level/verbosity
        self.V = DEBUG
        self.logger = logging.getLogger(__name__)
        if self.V == 1:
            lev = logging.DEBUG
        else:
            lev = logging.WARNING
        
        logging.basicConfig(format='%(asctime)s - [%(levelname)s]: %(message)s',
                            datefmt='%d-%b-%y %H:%M:%S', level=lev)
        
        self.logger.debug("initiating sigfish")
    
    def __init__(self, threads=8, DEBUG=0):
        '''
        python init
        '''
        somelist = []
    
    def __dealloc__(self):
        '''
        free memory
        '''
    

    def process_batch(self, batch):
        '''
        process a batch of of signals
        '''
        return