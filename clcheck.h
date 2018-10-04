// Remove all CL calls when compiling for CPU
#ifdef CL2CPU
	#define cl_int int
	#define cl_context int
	#define cl_device_id int
	#define cl_command_queue int
	#define cl_kernel int
	#define cl_mem int
	#define cl_event int
	#define cl_ulong unsigned long int
	#define CL(opname, ...)
	#define CLCREATE(opname, ...)

#else

	// Error reporting on each CL call when in DEBUG mode
	#ifdef DEBUG

		#include "Log.h"

		const char* clerr2str(cl_int);

		// internal macro
		#define CLCHECK_(opname) \
			if(clerrcode != CL_SUCCESS) { \
				out << clerr2str(clerrcode) << " returned by " << #opname; \
				out << " at line " << __LINE__ << " in file " << __FILE__ << '\n'; \
			}

		// user macros - prefix all your CL calls!
		#define CL(opname, ...) { cl_int clerrcode = cl ## opname (__VA_ARGS__); CLCHECK_(cl ## opname); }
		#define CLCREATE(opname, ...) { cl_int clerrcode = CL_SUCCESS; clCreate ## opname (__VA_ARGS__, &clerrcode); CLCHECK_(clCreate ## opname); }

	// Raw CL calls ignoring error codes when in non-debug mode
	#else

		// no checks in non-debug mode
		#define CL(opname, ...) cl ## opname (__VA_ARGS__)
		#define CLCREATE(opname, ...) clCreate ## opname (__VA_ARGS__, NULL)

	#endif

#endif