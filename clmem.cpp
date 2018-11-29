#include <assert.h>
#include <map>
#define DEBUG
#include "clcheck.h" // CL macro

void createGLBuffer(size_t size, void* outBuffer);

// Reusing CL memory handles:
// CL keeps the mapping to device memory internally
// Here we map the same handles to host pointers
static std::map<cl_mem, void*> map;

//TODO define a SAFE_MODE, where we keep a second map with array sizes, so we can crash the program with an error message before buffer overruns happen

// If no CL is present we use a simple counter as handle
#ifdef NO_GPU
static int handle = 0;
#endif

// Current OpenCL context that keeps the device memory pool
static cl_context context = 0;

// Buffers that have been shared with OpenGL
static cl_mem glObjects[10];
static int glObjectCount = 0;
int getCLMemSharedGLObjectCount() { return glObjectCount; }
cl_mem getCLMemSharedGLObject(int i) { return glObjects[i]; }

void setCLMemContext(cl_context c) {
	context = c;
}

cl_mem allocCLInputBuffer(size_t size) {
	#ifdef NO_GPU
	handle++;
	#else
	assert(context);
	cl_mem handle = CLCREATE(Buffer, context, CL_MEM_READ_ONLY, size, NULL);
	#endif
	void* pointer = malloc(size);
	assert(pointer);
	map[handle] = pointer;
	return handle;
}

cl_mem allocCLOutputBuffer(size_t size) {
	#ifdef NO_GPU
	handle++;
	#else
	assert(context);
	cl_mem handle;
	#ifdef GL_VISUALIZATION
	unsigned int glBuffer = 0;
	createGLBuffer(size, &glBuffer);
	if (glBuffer) {
		handle = CLCREATE(FromGLBuffer, context, CL_MEM_READ_WRITE, glBuffer);
		glObjects[glObjectCount++] = handle;
	} else handle = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, size, NULL);
	#else
	handle = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, size, NULL);
	#endif // GL_VISUALIZATION
	#endif // !NO_GPU
	void* pointer = malloc(size);
	assert(pointer);
	map[handle] = pointer;
	return handle;
}

void* getCLHostPointer(cl_mem handle) {
	return map[handle];
}

void freeCLMem() {
	for (const auto pair : map) {
		CL(ReleaseMemObject, pair.first);
		free(pair.second);
	}
}