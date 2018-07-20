#include <iostream> // std::cout
#include <string> // std::string
#include <stdint.h> // uint32_t
#include <stdlib.h> // malloc, free

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
		float nBelow = simulations[i].layers[simulations[i].n_layers + 1].n;
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
		int layerCount = simulations[i].n_layers;
		Layer* layers = (Layer*)malloc(layerCount * sizeof(Layer));
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

		int radialBinCount = simulations[i].det.nr;
		int angularBinCount = simulations[i].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint32_t);
		uint32_t* R_ra = (uint32_t*)malloc(reflectanceBufferSize);
		cl_mem gpuR_ra = CLCREATE(Buffer, context, CL_MEM_READ_WRITE, reflectanceBufferSize, NULL);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &gpuR_ra);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount);
		float radialBinCentimeters = simulations[i].det.dr;
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters);

		cl_event kernelEvent, reflectanceTransferEvent;
		CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
		CL(EnqueueReadBuffer, cmdQueue, gpuR_ra, CL_FALSE, 0, reflectanceBufferSize, R_ra, 0, NULL, &reflectanceTransferEvent);
		CL(Finish, cmdQueue);
		float totalDiffuseReflectance = 0.0f;
		for (int j = 0; j < radialBinCount; j++) {
			for (int k = 0; k < angularBinCount; k++) {
				float w = R_ra[j * angularBinCount + k] * (1.0f / 4294967296.0f);
				totalDiffuseReflectance += w;
			}
		}
		if (rank == 0) {
			cl_ulong timeStart, timeEnd;
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "PhotonCount=" << totalThreadCount << std::endl;
			std::cout << "Kerneltime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";
			std::cout << "totalDiffuseReflectance=" << totalDiffuseReflectance << std::endl;
		}
		free(R_ra);
		free(layers);
	}
	for (int i = 0; i < simCount; i++) {
		free(simulations[i].layers);
	}
	free(simulations);
}

// use struct of arrays for photon positions