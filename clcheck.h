#ifdef DEBUG

	static cl_int clerrcode = CL_SUCCESS;

	const char* clerr2str(cl_int errcode);

	#define CLCHECK(opname) \
		if(clerrcode != CL_SUCCESS) { \
			std::cout << clerr2str(clerrcode) << " returned by " << #opname; \
			std::cout << " at line " << __LINE__ << " in file " << __FILE__ << std::endl; \
		}

	#define CL(opname, ...) clerrcode = cl ## opname (__VA_ARGS__); CLCHECK(cl ## opname);

	#define CLCREATE(opname, ...) clCreate ## opname (__VA_ARGS__, &clerrcode); CLCHECK(clCreate ## opname);

#else
	// no checks in non-debug mode
		
	#define CL(opname, ...) cl ## opname (__VA_ARGS__)
	#define CLCREATE(opname, ...) clCreate ## opname (__VA_ARGS__, NULL)

#endif