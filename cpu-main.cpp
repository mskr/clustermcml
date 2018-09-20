#include <assert.h>

#define cl_int int
#define cl_context int
#define cl_command_queue int
#define cl_kernel int
#define cl_mem int
#define cl_event int
#define cl_ulong unsigned long int

void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
	int* inputBufferCount, size_t* inputBufferSizes,
	int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount);
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
	size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank);

int main(int nargs, char* args[]) {
	assert(nargs == 2);
	char* filename = args[1];
	{
		int inputBufferCount = 0, outputBufferCount = 0;
		size_t inputBufferSizes[10], outputBufferSizes[10];
		int inputBuffers[10], outputBuffers[10];
		allocCLKernelResources(1, "", filename, &inputBufferCount, inputBufferSizes, &outputBufferCount, outputBufferSizes, 10);
	}

	runCLKernel(0, 0, 0, 0, 0, 1, 1, 1, 0);
	return 0;
}