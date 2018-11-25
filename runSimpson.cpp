#include <iostream> // std::cout
#include <iomanip> // std::setprecision
#define DEBUG
#include "clcheck.h" // CL macro


static float* hostBuffer;


const char* getCLKernelName() {
	return "simpson";
}


void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
int* inputBufferCount, size_t* inputBufferSizes,
int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount, int rank) {

	// Report sizes for device buffers
	*inputBufferCount = 0;
	*outputBufferCount = 1;
	outputBufferSizes[0] = totalThreadCount * sizeof(float);

	// Alloc host buffers
	hostBuffer = new float[totalThreadCount]();
}


void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel simpsonKernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	
	// Global integration limits
	constexpr float a = 0;
	constexpr float b = 1;

	// Local integration limits
	float ha = a + (b - a) * rank/processCount;
	float hb = a + (b - a) * (rank+1)/processCount;

	cl_mem gpuBuffer = outputBuffers[0];

	CL(SetKernelArg, simpsonKernel, 0, sizeof(float), &ha);
	CL(SetKernelArg, simpsonKernel, 1, sizeof(float), &hb);
	CL(SetKernelArg, simpsonKernel, 2, sizeof(cl_mem), &gpuBuffer);

	cl_event simpsonKernelEvent;
	cl_event copyEvent;
	CL(EnqueueNDRangeKernel, cmdQueue, simpsonKernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &simpsonKernelEvent);
	CL(EnqueueReadBuffer, cmdQueue, gpuBuffer, CL_FALSE, 0, totalThreadCount * sizeof(float), hostBuffer, 0, NULL, &copyEvent);
	CL(Finish, cmdQueue);

	float value = 0;
	for (int i = 0; i < totalThreadCount; i++) {
		value += hostBuffer[i];
	}
	float sum;
	MPI_Reduce(&value, &sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (rank == 0) {
		cl_ulong timeSimpsonStart, timeSimpsonEnd;
		clGetEventProfilingInfo(simpsonKernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeSimpsonStart, NULL);
		clGetEventProfilingInfo(simpsonKernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeSimpsonEnd, NULL);
		std::cout << "Kernel time: " << (timeSimpsonEnd - timeSimpsonStart) << " ns\n";
		std::cout << "Used " << processCount << " processes and " << totalThreadCount << " threads for each." << std::endl;
		std::cout << std::setprecision(14) << sum << std::endl;
	}

	delete[] hostBuffer;
}