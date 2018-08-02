#include <iostream> // std::cout
#include <string> // std::string
#include <stdint.h> // uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

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

static SimulationStruct* simulations = 0;
static int simCount = 0;
static Layer** layersPerSimulation = 0;
static uint64_t** reflectancePerSimulation = 0;
static char* debugBuffer = 0;

void allocCLKernelResources(char* kernelOptions, char* mcmlOptions,
size_t* inputBufferCount, size_t* inputBufferSizes,
size_t* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount) {

	//TODO read with rank 0 process and broadcast input data
	int ignoreA = std::string(kernelOptions).find("-D IGNORE_A") != std::string::npos ? 1 : 0;
	std::cout << "--- "<<mcmlOptions<<" --->" << std::endl;
	simCount = read_simulation_data(mcmlOptions, &simulations, ignoreA);
	std::cout << "<--- "<<mcmlOptions<<" ---" << std::endl;

	assert(simCount <= maxBufferCount);

	layersPerSimulation = (Layer**)malloc(simCount * sizeof(Layer*));
	reflectancePerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));

	*inputBufferCount = simCount;
	*outputBufferCount = simCount;

	for (int i = 0; i < simCount; i++) {

		assert(simulations[i].number_of_photons <= 0xFFFFFFFFu); // ensures no bins can overflow

		int layerCount = simulations[i].n_layers;

		inputBufferSizes[i] = layerCount * sizeof(Layer);
		Layer* layers = (Layer*)malloc(inputBufferSizes[i]);
		layersPerSimulation[i] = layers;
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

		int radialBinCount = simulations[i].det.nr;
		int angularBinCount = simulations[i].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);
		outputBufferSizes[i] = reflectanceBufferSize;
		uint64_t* R_ra = (uint64_t*)malloc(reflectanceBufferSize);
		reflectancePerSimulation[i] = R_ra;
	}

	int debugMode = std::string(kernelOptions).find("-D DEBUG") != std::string::npos ? 1 : 0;
	if (debugMode) {
		debugBuffer = (char*)malloc(2048);
		outputBufferSizes[(*outputBufferCount)++] = 2048;
	}
}

static void freeResources() {
	for (int i = 0; i < simCount; i++) {
		free(reflectancePerSimulation[i]);
		free(layersPerSimulation[i]);
		free(simulations[i].layers);
	}
	free(reflectancePerSimulation);
	free(layersPerSimulation);
	free(simulations);
	if (debugBuffer) {
		free(debugBuffer);
	}
}

static bool handleDebugOutput() {
	const char* error = "error";
	bool isError = false;
	int j = 0;
	// print printable ascii chars if starting with "error"
	for (; debugBuffer[j] >= 32 && debugBuffer[j] <= 126; j++) {
		if (j <= 4 && debugBuffer[j] != error[j]) break;
		if (j == 4) {
			std::cout << error;
			isError = true;
		}
		if (j > 4) std::cout << debugBuffer[j];
	}
	if (isError) {
		std::cout << std::endl;
		for (int k = 0; k < 3; k++) { // print as floats
			std::cout << ((float*)debugBuffer)[j+k] << " ";
		}
		std::cout << std::endl;
	}
	return isError;
}

static void processReflectance(uint64_t* reflectanceArray, int radialBinCount, int angularBinCount,
double* Rd_r) {
	for (int i_r = 0; i_r < radialBinCount; i_r++) {
		uint64_t sumAngles = 0;
		for (int i_a = 0; i_a < angularBinCount; i_a++) {
			sumAngles += reflectanceArray[i_r * angularBinCount + i_a];
		}
		Rd_r[i_r] = (double)sumAngles / 0xFFFFFFFFf;
	}
}

//TODO use 64 bit bins
//TODO detect transmission and absorption
//TODO write output file
//TODO run original mcml
//TODO compare reflectance
//TODO simulate exact number of photons

void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {
		int argCount = 0;

		// first and last layer have only n initialized and represent outer media
		float nAbove = simulations[simIndex].layers[0].n;
		float nBelow = simulations[simIndex].layers[simulations[simIndex].n_layers + 1].n;

		int radialBinCount = simulations[simIndex].det.nr;
		float radialBinCentimeters = simulations[simIndex].det.dr;
		int angularBinCount = simulations[simIndex].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);

		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &inputBuffers[simIndex]);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &simulations[simIndex].n_layers);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[simIndex]);
		if (debugBuffer) {
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[simCount]);
		}

		cl_event kernelEvent, reflectanceTransferEvent;
		CL(EnqueueWriteBuffer, cmdQueue, inputBuffers[simIndex], CL_FALSE, 
			0, simulations[simIndex].n_layers * sizeof(Layer), layersPerSimulation[simIndex], 0, NULL, NULL);
		//TODO init accumulation buffers with zeros
		CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
		CL(EnqueueReadBuffer, cmdQueue, outputBuffers[simIndex], CL_FALSE, 
			0, reflectanceBufferSize, reflectancePerSimulation[simIndex], 0, NULL, &reflectanceTransferEvent);
		if (debugBuffer) {
			CL(EnqueueReadBuffer, cmdQueue, outputBuffers[simCount], CL_FALSE, 0, 2048, debugBuffer, 0, NULL, NULL);
		}
		CL(Finish, cmdQueue);

		if (debugBuffer) {
			if (handleDebugOutput()) return;
			// double sum = 0;
			// for (int i = 0; i < 512; i++) sum += (double)(((uint32_t*)debugBuffer)[i]) / 0xFFFFFFFF;
			// std::cout << (sum / 512.0) << std::endl;
		}

		if (rank == 0) {
			std::cout << "PhotonCount=" << totalThreadCount << std::endl;
			cl_ulong timeStart, timeEnd;
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "Kerneltime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";

			uint64_t* Rd_ra = reflectancePerSimulation[simIndex];

			Write_Simulation_Results(Rd_ra, Rd_ra, Rd_ra, &simulations[simIndex], timeEnd - timeStart);

			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "Transfertime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";
		}

		CL(ReleaseEvent, kernelEvent); CL(ReleaseEvent, reflectanceTransferEvent);
	}
	freeResources();
}