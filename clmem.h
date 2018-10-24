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

/**
* Data Problem 2
* --------------
* Handles instead of pointers for host memory.
* Safer through centralized memory management.
* Owner is always only this module.
* Handles also valid for corresponding device memory.
*/
#define CLMALLOC_INPUT(N, T) allocCLInputBuffer(N * sizeof(T));
#define CLMALLOC_OUTPUT(N, T) allocCLOutputBuffer(N * sizeof(T));
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

#endif // __OPENCL_VERSION__


/**
* Data Problem 3
* --------------
* Make memory layout explicit
* and easy to change.
* AOS == Array of Structs
* SOA == Struct of Arrays
*/

#if __OPENCL_VERSION__
#define CLMEM_ACCESS_ARRAY(pointer, i) (pointer[i])
#define CLMEM_ACCESS_ARRAY2D(pointer, size_j, i, j) (pointer[i * size_j + j])
#define CLMEM_ACCESS_AOS(pointer, i, member) (pointer[i].member)
#define CLMEM_ACCESS_SOA(pointer, i, member) (pointer->member[i])
#else
#define CLMEM_ACCESS_ARRAY(handle, T, i) (((T*)getCLHostPointer(handle))[i])
#define CLMEM_ACCESS_ARRAY2D(handle, T, size_j, i, j) (((T*)getCLHostPointer(handle))[i * size_j + j])
#define CLMEM_ACCESS_AOS(handle, T, i, member) (((T*)getCLHostPointer(handle))[i].member)
#define CLMEM_ACCESS_SOA(handle, T, i, member) (((T*)getCLHostPointer(handle))->member[i])
#endif