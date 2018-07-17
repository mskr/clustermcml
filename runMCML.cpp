#include <iostream> // std::cout
#include <string> // std::string

#define DEBUG
#include "clcheck.h"

#include "CUDAMCMLio.c"


struct Boundary {
	float z;
	float nx, ny, nz;
};

struct Layer {
	float absorbCoeff;
	float scatterCoeff;
	float g; // anisotropy
	float n; // refractive index
	struct Boundary top;
	struct Boundary bottom;
};

const char* getCLKernelName() {
	return "mcml";
}

void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank, char* kernelOptions, char* mcmlOptions) {
	int ignoreA = std::string(kernelOptions).find("-D IGNORE_A") != std::string::npos ? 1 : 0;
	SimulationStruct* simulations;
	std::cout << "--- start reading mcml input file "<<mcmlOptions<<" --->" << std::endl;
	int simCount = read_simulation_data(mcmlOptions, &simulations, ignoreA);
	std::cout << "<--- finished reading mcml input file "<<mcmlOptions<<" ---" << std::endl;
	for (int i = 0; i < simCount; i++) {
		int argCount = 0;
		// first and last layer have only n initialized and represent outer media
		float nAbove = simulations[i].layers[0].n;
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
		int layerCount = simulations[i].n_layers;
		float nBelow = simulations[i].layers[layerCount + 1].n;
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
		Layer* layers = new Layer[layerCount];
		// Q: why do mcml authors set scattering coeff which is probability to greater 1?
		for (int j = 1; j <= layerCount; j++) {
			layers[j - 1] = {
				simulations[i].layers[j].mua,
				1.0f / simulations[i].layers[j].mutr - simulations[i].layers[j].mua,
				simulations[i].layers[j].g,
				simulations[i].layers[j].n,
				Boundary{simulations[i].layers[j].z_min, 0.0f, 0.0f, 1.0f},
				Boundary{simulations[i].layers[j].z_max, 0.0f, 0.0f, -1.0f},
			};
		}
		cl_mem gpuLayers = CLCREATE(Buffer, context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, layerCount * sizeof(Layer), layers);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &gpuLayers);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &layerCount);
		unsigned int reflectCount = 0;
		unsigned int transmitCount = 0;
		unsigned int absorbCount = 0;
		cl_mem gpuReflectCount = CLCREATE(Buffer, context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &reflectCount);
		cl_mem gpuTransmitCount = CLCREATE(Buffer, context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &transmitCount);
		cl_mem gpuAbsorbCount = CLCREATE(Buffer, context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &absorbCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &gpuReflectCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &gpuTransmitCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &gpuAbsorbCount);
		cl_event event;
		CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &event);
		CL(EnqueueReadBuffer, cmdQueue, gpuReflectCount, CL_FALSE, 0, sizeof(unsigned int), &reflectCount, 0, NULL, NULL);
		CL(EnqueueReadBuffer, cmdQueue, gpuTransmitCount, CL_FALSE, 0, sizeof(unsigned int), &transmitCount, 0, NULL, NULL);
		CL(EnqueueReadBuffer, cmdQueue, gpuAbsorbCount, CL_FALSE, 0, sizeof(unsigned int), &absorbCount, 0, NULL, NULL);
		CL(Finish, cmdQueue);
		if (rank == 0) {
			cl_ulong timeStart, timeEnd;
			clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "PhotonCount=" << totalThreadCount << std::endl;
			// std::cout << "Kernel time: " << (timeEnd - timeStart) << " ns\n";
			std::cout << "ReflectCount=" << reflectCount << std::endl;
			std::cout << "TransmitCount=" << transmitCount << std::endl;
			std::cout << "AbsorbCount=" << absorbCount << std::endl;
		}
		delete[] layers;
	}
	for (int i = 0; i < simCount; i++) {
		free(simulations[i].layers);
	}
	free(simulations);
}

// use struct of arrays for photon positions