#ifdef NO_GPU
	// Currently in the NO_GPU build, also MPI is disabled,
	// since the only purpose right now is simpler debugging.
	// This can change in the future if MPI without CL is needed.
	#define MPI_Datatype int
	#define MPI_UINT32_T 0
	#define MPI_CHAR 0
	#define MPI_FLOAT 0
	#define MPICHECK(result)
	#define MPI(opname, ...)
#else
	#ifdef DEBUG

		#include "Log.h"

		#define MPICHECK(result) \
			if(result != MPI_SUCCESS) { \
				char err[64]; int len; \
				MPI_Error_string(result, err, &len); \
				out << err << " returned at line " << __LINE__ << " in file " << __FILE__ << '\n'; \
			}

		#define MPI(opname, ...) { int result = MPI_ ## opname (__VA_ARGS__); MPICHECK(result); }

	#else
		// no checks in non-debug mode
		#define MPICHECK(result) result
		#define MPI(opname, ...) MPI_ ## opname (__VA_ARGS__)
	#endif
#endif