/*********************************************************************************
*
* An assert macro for OpenCL kernels that writes error message to host buffer
* and then performs a devide-by-zero to crash the kernel.
*
*********************************************************************************/

//TODO dump the whole stack frame when assertions fail

#ifdef DEBUG
#define DEBUG_BUFFER debugBuffer
#define DEBUG_BUFFER_ARG , volatile __global char* DEBUG_BUFFER
#define STR_COPY(src, dst) for(int i=0; src[i]!='\0';i++) dst[i]=src[i];
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x) // The extra level of indirection 
// will allow the preprocessor to expand the macros 
// before they are converted to strings.
#define classert(expr)\
	if(!(expr)) {\
		const __constant char* _msg_ = "error: assertion failed at line " STR(__LINE__);\
		STR_COPY(_msg_, DEBUG_BUFFER);\
		DEBUG_BUFFER[0] = 42/0;\
	}
#else // no assert in non-debug mode
#define DEBUG_BUFFER_ARG
#define classert(expr)
#endif