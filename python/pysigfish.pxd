#cython: language_level=3
from libc.stdio cimport *
from libc.stdint cimport *
from libc.stdlib cimport *

cdef extern from "sigfish.h":

	cdef enum sigfish_status:
		SIGFISH_MORE = 0,      #more data needed
		SIGFISH_REJECT = 1,    #reject the read
		SIGFISH_CONT = 2        #continue with the read

	ctypedef struct sigfish_state_t:
		pass

	ctypedef struct sigfish_read_t:
		int32_t read_number;
		int32_t channel;
		uint64_t len_raw_signal;
		float* raw_signal;
		pass
	
	# sigfish interface
	sigfish_state_t *init_sigfish(const char *ref, int num_channels, int threads);
	sigfish_status *process_sigfish(sigfish_state_t *state, sigfish_read_t *read_batch, int batch_size);
	void free_sigfish(sigfish_state_t *state);