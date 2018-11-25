#include <assert.h> // assert

#define DEBUG
#include "clcheck.h"

#ifndef CL2CPU
#include "clusterlib.h"
#endif

// interfaces for external code
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
	int* inputBufferCount, size_t* inputBufferSizes,
	int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount, int rank);
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
	size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank);
#ifdef GL_VISUALIZATION
void createGLBuffer(size_t size, void* outBuffer);
void runGLRenderLoop();
#endif

/**
*
*/
int main(int nargs, char** args) {

	#ifdef CL2CPU

	// Expect input file name
	assert(nargs == 2);
	char* filename = args[1];

	// Bypass all multiprocessor stuff and use dummy info
	int inputBufferCount = 0, outputBufferCount = 0;
	size_t inputBufferSizes[10], outputBufferSizes[10];
	int inputBuffers[10], outputBuffers[10];
	allocCLKernelResources(1, "", filename, &inputBufferCount, inputBufferSizes, &outputBufferCount, outputBufferSizes, 10, 0);
	runCLKernel(0, 0, 0, 0, 0, 1, 1, 1, 0);
	return 0;

	#else

	// Init and get cluster info (processes, threads etc.)
	int processCount;
	cl_context context;
	cl_command_queue cmdQueue1;
	cl_kernel kernel;
	size_t totalThreadCount, simdThreadCount;
	initCluster(nargs, args, &processCount, &context, &cmdQueue1, &kernel, &totalThreadCount, &simdThreadCount);

	// Get kernel configuration (defines etc.)
	char* kernelOptions = nargs >= 3 ? args[2] : "";

	// Get application specific command line options (input filenames etc.)
	char* otherOptions = nargs >= 4 ? args[3] : "";

	// Get application specific memory requirements
	int inputBufferCount = 0, outputBufferCount = 0;
	size_t inputBufferSizes[10], outputBufferSizes[10];
	cl_mem inputBuffers[10], outputBuffers[10];
	allocCLKernelResources(totalThreadCount, kernelOptions, otherOptions, &inputBufferCount, inputBufferSizes, &outputBufferCount, outputBufferSizes, 10, getRank());

	// Host buffers were created and sizes reported, now create CL buffers, optionally with GL interop
	for (int i = 0; i < inputBufferCount; i++) {
		inputBuffers[i] = CLCREATE(Buffer, context, CL_MEM_READ_ONLY, inputBufferSizes[i], NULL);
	}
	#ifdef GL_VISUALIZATION
	bool* couldCreateGLBuffer = new bool[outputBufferCount]();
	#endif
	for (int i = 0; i < outputBufferCount; i++) {
		#ifdef GL_VISUALIZATION
		unsigned int glBuf = 0;
		createGLBuffer(outputBufferSizes[i], &glBuf);
		if (glBuf) {
			couldCreateGLBuffer[i] = true;
			outputBuffers[i] = CLCREATE(FromGLBuffer, context, CL_MEM_READ_WRITE, glBuf);
			CL(EnqueueAcquireGLObjects, cmdQueue1, 1, &outputBuffers[i], 0, NULL, NULL);
		} else {
			outputBuffers[i] = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, outputBufferSizes[i], NULL);
		}
		#else
		outputBuffers[i] = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, outputBufferSizes[i], NULL);
		#endif
	}

	// Run actual application
	runCLKernel(context, cmdQueue1, kernel, inputBuffers, outputBuffers, totalThreadCount, simdThreadCount, processCount, getRank());

	// Release objects so that they can be used by OpenGL
	#ifdef GL_VISUALIZATION
	for (int i = 0; i < outputBufferCount; i++) {
		if (couldCreateGLBuffer[i]) {
			CL(EnqueueReleaseGLObjects, cmdQueue1, 1, &outputBuffers[i], 0, NULL, NULL);
		}
	}
	delete[] couldCreateGLBuffer;
	runGLRenderLoop();
	#endif

	cleanupCluster(inputBufferCount, inputBuffers, outputBufferCount, outputBuffers, context, cmdQueue1, kernel);

	#endif // CL2CPU
}

// hello world example
// http://phycomp.technion.ac.il/~comphy/classfiles/hellocl.c

// nvidia optimization: occupancy, ...
// https://www.cs.cmu.edu/afs/cs/academic/class/15668-s11/www/cuda-doc/OpenCL_Best_Practices_Guide.pdf

// Single kernel is for data-parallel problems, running multiple kernels in parallel is for task-parallel problems
// Can use multiple queues or out-of-order-mode to enqueue multiple kernels
// https://us.fixstars.com/opencl/book/OpenCLProgrammingBook/opencl-programming-practice/

// work can be further spread on multiple devices using multiple command queues

// For also using CPU there are options:
// 1. create a separate cl context and compile kernel for CPU
// 2. compile kernel with C compiler (macro defines, no vector types) 
// 3. Automatic scheduling of the submitted commands on OpenCL devices
//    http://starpu.gforge.inria.fr/doc/html/SOCLOpenclExtensions.html

// Hot reloading hardcoded parameters
// http://tuxedolabs.blogspot.com/2018/03/hot-reloading-hardcoded-parameters.html?m=1
// File watcher: https://facebook.github.io/watchman/