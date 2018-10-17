/**
* Data Problem 1
* --------------
* Explicit and cross-platform alignment of structs
* ALIGN_X(
* struct A {
*   int a;
* });
*/
#if _MSC_VER && !__INTEL_COMPILER // MSVC
#define ALIGN_NONE(T) __pragma(pack(push, 1)) T __pragma(pack(pop))
#define ALIGN_4BYTE(T) __pragma(pack(push, 4)) T __pragma(pack(pop))
#define ALIGN_8BYTE(T) __pragma(pack(push, 8)) T __pragma(pack(pop))
#else // GCC or OpenCL
#define ALIGN_NONE(T) T __attribute__((aligned(1)))
#define ALIGN_4BYTE(T) T __attribute__((aligned(4)))
#define ALIGN_8BYTE(T) T __attribute__((aligned(8)))
#endif

#define CLMALLOC(N, T) allocCLMem(N * sizeof(T));

/**
* Data Problem 2
* --------------
* Handles instead of pointers for host memory.
* Safer and also usable as handles
* for corresponding device memory.
*/
cl_mem allocCLInputBuffer(size_t size);
cl_mem allocCLOutputBuffer(size_t size);

void* getCLMemHostPointer(cl_mem handle);


/**
* Data Problem 3
* --------------
* Easily change memory layouts
* via accessor functions.
*/

// For layout specific accessors we either
// need a language feature to inspect structs or
// write a method for each datatype or
// use templates which we cannot compile for GPU.
// The problem of switching data representations
// by changing code only in one place
// seems impossible to solve in C.

enum Layout { ARRAY, ARRAY2D, AOS, SOA };

void setCLMemAsAOS(void* pointer, size_t index, Layer data);