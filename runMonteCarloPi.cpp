#include <iostream> // std::cout
#include <iomanip> // std::setprecision

const char* getCLKernelName() {
	return "mcpi";
}

void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	unsigned int* hostBuffer = new unsigned int[totalThreadCount];
	cl_mem gpuBuffer = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, totalThreadCount * sizeof(unsigned int), NULL);
	int npoints = 1000;
	CL(SetKernelArg, kernel, 0, sizeof(int), &npoints);
	CL(SetKernelArg, kernel, 1, sizeof(cl_mem), &gpuBuffer);
	cl_event kernelEvent;
	CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
	CL(EnqueueReadBuffer, cmdQueue, gpuBuffer, CL_FALSE, 0, totalThreadCount * sizeof(unsigned int), hostBuffer, 0, NULL, NULL);
	CL(Finish, cmdQueue);
	unsigned long sum = 0;
	for (int i = 0; i < totalThreadCount; i++) {
		sum += hostBuffer[i];
	}
	unsigned long totalSum = 0;
	MPI_Reduce(&sum, &totalSum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	// 1. have a unit square (area 1) with a quarter circle (area pi/4) inside
	// 2. shoot n random points uniformly at the square
	// 3. when h points hit the circle => (pi/4)/1 = h/n <=> pi = 4*h/n
	long totalPoints = npoints * totalThreadCount * processCount;
	double pi = 4.0 * totalSum / totalPoints;
	if (rank == 0) {
		cl_ulong timeStart, timeEnd;
		clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
		clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
		std::cout << "Kernel time: " << (timeEnd - timeStart) << " ns\n";
		std::cout << "totalPoints=" << totalPoints << ", hitPoints=" << totalSum << std::endl;
		std::cout << std::setprecision(14) << pi << std::endl;
	}
}