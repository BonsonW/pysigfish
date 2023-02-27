# distutils: language = c
# cython: language_level=3
# cython: profile=True
import sys
import time
import logging
import copy
from libc.stdlib cimport malloc, free
from libc.string cimport strdup
cimport pysigfish
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
    cdef sigfish_state_t *state
    cdef char* REF
    cdef int NUM_CHANNELS
    cdef int NUM_THREADS
    cdef sigfish_read_t *sbatch
    cdef sigfish_status *status
    cdef int batch_len


    def __cinit__(self, ref, channels=512, threads=8, DEBUG=0):
        '''
        C init
        '''
        self.state = NULL
        self.REF = ""
        self.NUM_CHANNELS = 0
        self.NUM_THREADS = 0
        self.sbatch = NULL
        self.status = NULL
        self.batch_len = 0


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

        REF = str.encode(ref)
        self.REF = strdup(REF)
        self.NUM_CHANNELS = channels
        self.NUM_THREADS = threads

        self.state = init_sigfish(self.REF, self.NUM_CHANNELS, self.NUM_THREADS)
        if self.state is NULL:
            self.logger.error("Ref '{}' could not be opened and sigfish not initialised".format(ref))
        
        if self.state is NULL:
            raise MemoryError()


    
    def __init__(self, ref, channels=512, threads=8, DEBUG=0):
        '''
        python init
        '''
        somelist = []
    
    def __dealloc__(self):
        '''
        free memory
        '''
        if self.state is not NULL:
            free_sigfish(self.state)

        if self.REF is not NULL:
            free(self.REF)
    

    def process_batch(self, batch, signal_dtype):
        '''
        process a batch of of signals
        ctypedef struct sigfish_read_t:
            int32_t read_number;
            int32_t channel;
            uint64_t len_raw_signal;
            float* raw_signal;
        
        cdef enum sigfish_status:
            SIGFISH_MORE = 0,      #more data needed
            SIGFISH_REJECT = 1,    #reject the read
            SIGFISH_CONT = 2        #continue with the read
        
        readID = read.id
        read_number = read.number
        chunk_length = read.chunk_length
        raw_data = numpy.fromstring(read.raw_data, dtype)
        '''
        self.batch_len = len(batch)
        self.sbatch = <sigfish_read_t *> malloc(sizeof(sigfish_read_t)*self.batch_len)

        for idx, channel, read in enumerate(batch):
            self.sbatch[idx].read_number = read.number
            self.sbatch[idx].channel = channel
            self.sbatch[idx].len_raw_signal = read.chunk_length
            self.sbatch[idx].raw_signal = <float *> malloc(sizeof(float)*read.chunk_length)
            sig = np.fromstring(read.raw_data, signal_dtype)
            memview = memoryview(sig)
            for i in range(read.chunk_length):
                self.sbatch[idx].raw_signal[i] = memview[i]

        self.status = process_sigfish(self.state, self.sbatch, self.batch_len)

        status_dic = {}
        for idx, channel, read in enumerate(batch):
            status_dic[channel] = (read.channel, read.number, read.id, self.status[idx])
        
        # free memory
        for i in range(self.batch_len):
            for j in range(self.sbatch[i].len_raw_signal):
                free(self.sbatch[i].raw_signal)
        free(self.sbatch)
        free(self.status)

        return status_dic

