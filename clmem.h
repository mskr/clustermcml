/*********************************************************************************
*
* This file serves as centralized memory manager supporting CPU and GPU code.
*
* Centralized memory management is considered more safe than using pointers or
* even smart pointers all over the code.
*
* In detail, it takes care of
* 1) aligning struct data in a cross-platform way for easy interchange,
* 2) allocating host and device memory in one go returning a unified handle,
* 3) accessing memory with explicit declaration of memory layout (including easy 
*    switch between Array of Structs (AOS) to Struct of Arrays (SOA)).
*
*********************************************************************************/

/**
* Data Problem 1
* --------------
* Explicit and cross-platform alignment of structs
* ALIGN_X(
* struct A {
*   int a;
* });
*/
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER) // MSVC
#define ALIGN_NONE(T) __pragma(pack(push, 1)) T __pragma(pack(pop))
#define ALIGN_4BYTE(T) __pragma(pack(push, 4)) T __pragma(pack(pop))
#define ALIGN_8BYTE(T) __pragma(pack(push, 8)) T __pragma(pack(pop))
#elif defined(__GNUC__) || defined(__OPENCL_VERSION__) // GCC or OpenCL
#define ALIGN_NONE(T) T __attribute__((aligned(1)))
#define ALIGN_4BYTE(T) T __attribute__((aligned(4)))
#define ALIGN_8BYTE(T) T __attribute__((aligned(8)))
#endif


#if !defined(__OPENCL_VERSION__)

#ifdef NO_GPU
#define cl_mem int
#define cl_context int
#endif

/**
* Data Problem 2
* --------------
* Handles instead of pointers for host memory.
* Safer through centralized memory management.
* Owner is always only this module.
* Handles also valid for corresponding device memory.
*/
#define CLMALLOC_INPUT(N, T) allocCLInputBuffer(N * sizeof(T))
#define CLMALLOC_OUTPUT(N, T) allocCLOutputBuffer(N * sizeof(T))
#define CLMEM(handle) getCLHostPointer(handle)
cl_mem allocCLInputBuffer(size_t size);
cl_mem allocCLOutputBuffer(size_t size);
void* getCLHostPointer(cl_mem handle);

/**
* Set OpenCL context for which allocations will happen.
* OpenCL shares memory for multiple devices via context.
*/
void setCLMemContext(cl_context);

/**
* Get output buffers that have been shared with OpenGL.
* Sharing is done when compiled with GL_VISUALIZATION define.
*/
int getCLMemSharedGLObjectCount();
cl_mem getCLMemSharedGLObject(int i);

void freeCLMem();

#endif // !defined(__OPENCL_VERSION__)


/**
* Data Problem 3
* --------------
* Make memory layout explicit
* and easy to change.
* AOS == Array of Structs
* SOA == Struct of Arrays
*/

#define CLMEM_ACCESS_ARRAY(pointer, T, i) (((T*)pointer)[i])
#define CLMEM_ACCESS_ARRAY2D(pointer, T, size_j, i, j) (((T*)pointer)[i * size_j + j])
#define CLMEM_ACCESS_AOS(pointer, T, i, member) (((T*)pointer)[i].member)
#define CLMEM_ACCESS_SOA(pointer, T, i, member) (((T*)pointer)->member[i])