#include <iostream> // std::cout
#include <iomanip> // std::setprecision
#define DEBUG
#include "clcheck.h" // CL macro

static unsigned int* hostBuffer;


const char* getCLKernelName() {
	return "mcpi";
}


void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
int* inputBufferCount, size_t* inputBufferSizes,
int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount) {

	// Report sizes for device buffers
	*inputBufferCount = 0;
	*outputBufferCount = 1;
	outputBufferSizes[0] = totalThreadCount * sizeof(unsigned int);

	// Alloc host buffers
	hostBuffer = new unsigned int[totalThreadCount];
}


void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {

	// 1. have a unit square (area 1) with a quarter circle (area pi/4) inside
	// 2. shoot n random points uniformly at the square
	// 3. when h points hit the circle => (pi/4)/1 = h/n <=> pi = 4*h/n

	// Number of random points
	int npoints = 1000;

	// Get device buffer
	cl_mem gpuBuffer = outputBuffers[0];

	// Set args
	CL(SetKernelArg, kernel, 0, sizeof(int), &npoints);
	CL(SetKernelArg, kernel, 1, sizeof(cl_mem), &gpuBuffer);

	// Run kernel
	cl_event kernelEvent;
	CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);

	// Download device buffer
	CL(EnqueueReadBuffer, cmdQueue, gpuBuffer, CL_FALSE, 0, totalThreadCount * sizeof(unsigned int), hostBuffer, 0, NULL, NULL);
	CL(Finish, cmdQueue);

	// Collect results from device buffer
	unsigned long sum = 0;
	for (int i = 0; i < totalThreadCount; i++) {
		sum += hostBuffer[i];
	}

	// Collect results from all processes
	unsigned long totalSum = 0;
	MPI_Reduce(&sum, &totalSum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	// Calculate PI
	long totalPoints = npoints * totalThreadCount * processCount;
	double pi = 4.0 * totalSum / totalPoints;

	// Write output
	if (rank == 0) {
		cl_ulong timeStart, timeEnd;
		clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
		clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
		std::cout << "Kernel time: " << (timeEnd - timeStart) << " ns\n";
		std::cout << "totalPoints=" << totalPoints << ", hitPoints=" << totalSum << std::endl;
		std::cout << std::setprecision(14) << pi << std::endl;
	}
}