#ifdef NO_GPU
#include <stdint.h>
typedef uint32_t uint;
#define __global
#elif !defined(__OPENCL_VERSION__)
typedef cl_float3 float3; // note: sizeof(cl_float3) != sizeof(float[3])
typedef cl_uint uint;
#define __global
#endif
#include "geometrylib.h" // struct RHeightfield

ALIGN_4BYTE(
struct Boundary {
	
	uint isHeightfield;

	// valid if !isHeightfield
	float z;

	// Padding that is implicitely added by OpenCL compiler,
	// but not by regular C++ compiler
	//TODO solve this problem with alignment attribute or member ordering
	float padding1;
	float padding2;

	// valid if isHeightfield
	struct RHeightfield heightfield;
});