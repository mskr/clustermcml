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
	int simCount = read_simulation_data(mcmlOptions, &simulations, ignoreA);
	for (int i = 0; i < simCount; i++) {
		float outerN = 1.000292f; // air
		CL(SetKernelArg, kernel, 0, sizeof(float), &outerN);
		int layerCount = simulations[i].n_layers;
		Layer* layers = new Layer[layerCount];
		// Q: why do mcml authors set scattering coeff which is probability to greater 1?
		for (int j = 0; j < layerCount; j++) {
			layers[j] = {
				simulations[i].layers[j].mua,
				1.0f / simulations[i].layers[j].mutr - simulations[i].layers[j].mua,
				simulations[i].layers[j].g,
				simulations[i].layers[j].n,
				Boundary{simulations[i].layers[j].z_min, 0.0f, 0.0f, 1.0f},
				Boundary{simulations[i].layers[j].z_max, 0.0f, 0.0f, -1.0f},
			};
			std::cout << simulations[i].layers[j].z_min << " " << simulations[i].layers[j].z_max << std::endl; //TODO read unreasonable z values
		}
		cl_mem gpuLayers = CLCREATE(Buffer, context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, layerCount * sizeof(Layer), layers);
		CL(SetKernelArg, kernel, 1, sizeof(cl_mem), &gpuLayers);
		CL(SetKernelArg, kernel, 2, sizeof(int), &layerCount);
		unsigned int reflectCount = 0;
		unsigned int transmitCount = 0;
		unsigned int absorbCount = 0;
		cl_mem gpuReflectCount = CLCREATE(Buffer, context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &reflectCount);
		cl_mem gpuTransmitCount = CLCREATE(Buffer, context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &transmitCount);
		cl_mem gpuAbsorbCount = CLCREATE(Buffer, context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int), &absorbCount);
		CL(SetKernelArg, kernel, 3, sizeof(cl_mem), &gpuReflectCount);
		CL(SetKernelArg, kernel, 4, sizeof(cl_mem), &gpuTransmitCount);
		CL(SetKernelArg, kernel, 5, sizeof(cl_mem), &gpuAbsorbCount);
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