#ifdef DEBUG
	// have these global variables only in debug mode
	static int rank = 0;
	static cl_int clerrcode = CL_SUCCESS;
	const char* clerr2str(cl_int errcode);
	#define MPICHECK(result) \
		if(result != MPI_SUCCESS) { \
			char err[64]; int len; \
			MPI_Error_string(result, err, &len); \
			std::cout << err << " returned at line " << __LINE__ << " in file " << __FILE__ << std::endl; \
		}
	#define CLLOG(opname) \
		if(clerrcode != CL_SUCCESS && rank == 0) { \
			std::cout << clerr2str(clerrcode) << " returned by " << #opname; \
			std::cout << " at line " << __LINE__ << " in file " << __FILE__ << std::endl; \
		}
	#define CL(opname, ...) clerrcode = cl ## opname (__VA_ARGS__); CLLOG(cl ## opname);
	#define CLCREATE(opname, ...) clCreate ## opname (__VA_ARGS__, &clerrcode); CLLOG(clCreate ## opname);
#else
	#define MPICHECK(result) result
	#define CL(opname, ...) cl ## opname (__VA_ARGS__)
	#define CLCREATE(opname, ...) clCreate ## opname (__VA_ARGS__, NULL)
#endif