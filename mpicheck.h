#ifdef DEBUG

	#define MPICHECK(result) \
		if(result != MPI_SUCCESS) { \
			char err[64]; int len; \
			MPI_Error_string(result, err, &len); \
			std::cout << err << " returned at line " << __LINE__ << " in file " << __FILE__ << std::endl; \
		}

#else
	// no checks in non-debug mode
	#define MPICHECK(result) result
#endif