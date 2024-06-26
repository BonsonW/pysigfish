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
		char* read_id;
		pass
	
	ctypedef struct sigfish_opt_t:
		int num_thread;
		const char* debug_paf;
		int8_t no_full_ref;
		float dtw_cutoff;
		int query_size_sig;
		int query_size_events;
		int8_t pore;
	
	# sigfish interface
	sigfish_state_t *init_sigfish(const char *ref, int num_channels, sigfish_opt_t opt);
	sigfish_status *process_sigfish(sigfish_state_t *state, sigfish_read_t *read_batch, int batch_size);
	void free_sigfish(sigfish_state_t *state);

cdef extern from "error.h":

	# the level of verbosity in the log printed to the standard error
	cdef enum sigfish_log_level_opt:
		LOG_OFF=0,      # nothing at all
		LOG_ERR,        # error messages
		LOG_WARN,       # warning and error messages
		LOG_INFO,       # information, warning and error messages
		LOG_VERB,       # verbose, information, warning and error messages
		LOG_DBUG,       # debugging, verbose, information, warning and error messages
		LOG_TRAC        # tracing, debugging, verbose, information, warning and error messages
	
	void set_log_level(sigfish_log_level_opt level);
