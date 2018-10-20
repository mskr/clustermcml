#include <assert.h> // assert

#define DEBUG
#include "clcheck.h"

#ifndef NO_GPU
#include "clusterlib.h"
#endif

#include "clmem.h"

// interfaces for external code
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions, int rank);
void runCLKernel(cl_command_queue cmdQueue, cl_kernel kernel, size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank);
#ifdef GL_VISUALIZATION
void runGLRenderLoop();
#endif

/**
*
*/
int main(int nargs, char** args) {

	#ifdef NO_GPU

	// Expect input file name
	assert(nargs == 2);
	char* filename = args[1];

	// Bypass all multiprocessor stuff and use dummy info
	allocCLKernelResources(1, "", filename, 0);
	runCLKernel(0, 0, 1, 1, 1, 0);
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

	setCLMemContext(context);

	// Allocate resources based on application specific memory requirements
	allocCLKernelResources(totalThreadCount, kernelOptions, otherOptions, getRank());

	#ifdef GL_VISUALIZATION
	// Acquire shared GL buffers for use by CL from the device queue that does the visualization
	for (int i = 0; i < getCLMemSharedGLObjectCount(); i++)
		CL(EnqueueAcquireGLObjects, cmdQueue1, 1, &getCLMemSharedGLObject(i), 0, NULL, NULL);
	#endif

	// Run actual application
	runCLKernel(cmdQueue1, kernel, totalThreadCount, simdThreadCount, processCount, getRank());

	#ifdef GL_VISUALIZATION
	// Release shared GL buffers as they can now be used by GL
	for (int i = 0; i < getCLMemSharedGLObjectCount(); i++)
		CL(EnqueueReleaseGLObjects, cmdQueue1, 1, &getCLMemSharedGLObject(i), 0, NULL, NULL);
	}
	runGLRenderLoop();
	#endif

	freeCLMem();
	cleanupCluster(context, cmdQueue1, kernel);

	#endif // NO_GPU
}

// hello world example
// http://phycomp.technion.ac.il/~comphy/classfiles/hellocl.c

// Nvidia IHV Talk - Hardware and Optimizations (Timo Stich, Nvidia)
// http://sa09.idav.ucdavis.edu/docs/SA09_NVIDIA_IHV_talk.pdf

// AMD IHV Talk - Hardware and Optimizations (Jason Yang, AMD)
// http://sa09.idav.ucdavis.edu/docs/SA09_AMD_IHV.pdf

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
// 4. Maybe more in: Advanced OpenCL Event Model Usage by Derek Gerstmann
//    http://sa09.idav.ucdavis.edu/docs/SA09-opencl-dg-events-stream.pdf

// Hot reloading hardcoded parameters
// http://tuxedolabs.blogspot.com/2018/03/hot-reloading-hardcoded-parameters.html?m=1
// File watcher: https://facebook.github.io/watchman/

// Wave intrinsics from Microsoft shader compiler point of view
// https://github.com/Microsoft/DirectXShaderCompiler/wiki/Wave-Intrinsics