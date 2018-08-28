#include <iostream> // std::cout
#include <string> // std::string
#include <stdint.h> // uint32_t, uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

#define DEBUG
#include "clcheck.h"


//TODO share structs with kernel via header

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

struct PhotonState {
	float x, y, z; // pos [cm]
	float dx, dy, dz; // dir
	float weight; // 1 at start, zero when terminated
	int layerIndex; // current layer
};


static PhotonState createNewPhotonState() {
	return {
		0.0f, 0.0f, 0.0f, // start position
		0.0f, 0.0f, 1.0f, // start direction
		1.0f, // start weight
		0 // start layer index
	};
}


const char* getCLKernelName() {
	return "mcml";
}


static SimulationStruct* simulations = 0;
static int simCount = 0;
static Layer** layersPerSimulation = 0;
static uint64_t** reflectancePerSimulation = 0;
static uint64_t** transmissionPerSimulation = 0;
static uint64_t** absorptionPerSimulation = 0;
static PhotonState** photonStatesPerSimulation = 0;
static char* debugBuffer = 0;


void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* mcmlOptions,
int* inputBufferCount, size_t* inputBufferSizes,
int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount) {

	//TODO read with rank 0 and broadcast input data
	int ignoreA = std::string(kernelOptions).find("-D IGNORE_A") != std::string::npos ? 1 : 0;
	std::cout << "--- "<<mcmlOptions<<" --->" << std::endl;
	simCount = read_simulation_data(mcmlOptions, &simulations, ignoreA);
	std::cout << "<--- "<<mcmlOptions<<" ---" << std::endl;

	layersPerSimulation = (Layer**)malloc(simCount * sizeof(Layer*));
	reflectancePerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));
	transmissionPerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));
	absorptionPerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));
	photonStatesPerSimulation = (PhotonState**)malloc(simCount * sizeof(PhotonState*));

	*inputBufferCount = simCount;
	*outputBufferCount = simCount * 2; //TODO add T and A buffers

	assert(*outputBufferCount <= maxBufferCount);

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
		uint64_t* R_ra = (uint64_t*)malloc(reflectanceBufferSize);
		reflectancePerSimulation[i] = R_ra;

		size_t transmissionBufferSize = reflectanceBufferSize;
		uint64_t* T_ra = (uint64_t*)malloc(transmissionBufferSize);
		transmissionPerSimulation[i] = T_ra;

		int depthBinCount = simulations[i].det.nz;
		size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(uint64_t);
		uint64_t* A_rz = (uint64_t*)malloc(absorptionBufferSize);
		absorptionPerSimulation[i] = A_rz;

		photonStatesPerSimulation[i] = (PhotonState*)malloc(totalThreadCount * sizeof(PhotonState));

		outputBufferSizes[i] = reflectanceBufferSize;
		outputBufferSizes[simCount + i] = totalThreadCount * sizeof(PhotonState);
	}

	int debugMode = std::string(kernelOptions).find("-D DEBUG") != std::string::npos ? 1 : 0;
	if (debugMode) {
		debugBuffer = (char*)malloc(2048);
		outputBufferSizes[(*outputBufferCount)++] = 2048;
		assert(*outputBufferCount <= maxBufferCount);
	}
}


static void freeResources() {
	if (debugBuffer) {
		free(debugBuffer);
	}
	for (int i = 0; i < simCount; i++) {
		free(photonStatesPerSimulation[i]);
		free(reflectancePerSimulation[i]);
		free(layersPerSimulation[i]);
		free(simulations[i].layers);
	}
	free(reflectancePerSimulation);
	free(layersPerSimulation);
	free(simulations);
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
		//TODO print everything as hex view until reaching some unique end symbol
		for (int k = 0; k < 3; k++) { // print as floats
			std::cout << ((float*)debugBuffer)[j+k] << " ";
		}
		std::cout << std::endl;
	}
	return isError;
}


//TODO detect transmission and absorption
//TODO write output file
//TODO run original mcml
//TODO compare


void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		// first and last layer have only n initialized and represent outer media
		float nAbove = simulations[simIndex].layers[0].n;
		float nBelow = simulations[simIndex].layers[simulations[simIndex].n_layers + 1].n;

		// Reflectance (specular):
		// percentage of light leaving at surface without any interaction
		// using Fesnel approximation by Schlick (no incident angle, no polarization)
		// Q: why are the ^2 different than in Schlicks approximation?
		float nDiff = nAbove - simulations[simIndex].layers[0].n;
		float nSum = nAbove + simulations[simIndex].layers[0].n;
		float R_specular = (nDiff * nDiff) / (nSum * nSum);
		
		// init photon states
		for (int i = 0; i < totalThreadCount; i++) {
			PhotonState newState = createNewPhotonState();
			newState.weight -= R_specular;
			photonStatesPerSimulation[simIndex][i] = newState;
		}

		int radialBinCount = simulations[simIndex].det.nr;
		float radialBinCentimeters = simulations[simIndex].det.dr;
		int angularBinCount = simulations[simIndex].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);

		size_t photonStateBufferSize = totalThreadCount * sizeof(PhotonState);

		// init accumulation buffers with zeros
		for (int i = 0; i < radialBinCount * angularBinCount; i++) {
			reflectancePerSimulation[simIndex][i] = 0;
		}

		int argCount = 0;
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &inputBuffers[simIndex]);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &simulations[simIndex].n_layers);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount);
		CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters);
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[simIndex]); // reflectance buffer
		CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[simCount + simIndex]); // photon state buffer
		if (debugBuffer) {
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[simCount]);
		}

		cl_event kernelEvent, reflectanceTransferEvent;

		CL(EnqueueWriteBuffer, cmdQueue, inputBuffers[simIndex], CL_FALSE, 0,
			simulations[simIndex].n_layers * sizeof(Layer), layersPerSimulation[simIndex], 0, NULL, NULL);
		CL(EnqueueWriteBuffer, cmdQueue, outputBuffers[simIndex], CL_FALSE, 0,
			reflectanceBufferSize, reflectancePerSimulation[simIndex], 0, NULL, NULL);

		// Run kernel with optimal thread count as long as targeted number of photons allows it
		uint32_t targetPhotonCount = simulations[simIndex].number_of_photons;
		uint32_t finishedPhotonCount = 0;
		while (finishedPhotonCount < targetPhotonCount) { // stop when target reached
			CL(EnqueueWriteBuffer, cmdQueue, outputBuffers[simCount + simIndex], CL_FALSE, 0,
				photonStateBufferSize, photonStatesPerSimulation[simIndex], 0, NULL, NULL);
			size_t remainingPhotonCount = targetPhotonCount - finishedPhotonCount;
			if (remainingPhotonCount > totalThreadCount) {
				CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
			} else {
				CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &remainingPhotonCount, NULL, 0, NULL, &kernelEvent);
			}
			CL(EnqueueReadBuffer, cmdQueue, outputBuffers[simCount + simIndex], CL_FALSE, 0,
				photonStateBufferSize, photonStatesPerSimulation[simIndex], 0, NULL, NULL);
			if (debugBuffer) {
				CL(EnqueueReadBuffer, cmdQueue, outputBuffers[simCount], CL_FALSE, 0, 2048, debugBuffer, 0, NULL, NULL);
			}
			CL(Finish, cmdQueue);
			if (debugBuffer) {
				assert(!handleDebugOutput());
			}
			for (int i = 0; i < totalThreadCount; i++) {
				if (photonStatesPerSimulation[simIndex][i].weight == 0) {
					finishedPhotonCount++;
					PhotonState newState = createNewPhotonState();
					newState.weight -= R_specular;
					photonStatesPerSimulation[simIndex][i] = newState;
					//TODO since buffer updates are sparse, map could be faster than write in whole
				}
			}
			if (rank == 0) {
				std::cout << "Photons terminated: " << finishedPhotonCount << "/" << targetPhotonCount << std::endl;
			}
		}

		CL(EnqueueReadBuffer, cmdQueue, outputBuffers[simIndex], CL_FALSE, 0,
			reflectanceBufferSize, reflectancePerSimulation[simIndex], 0, NULL, &reflectanceTransferEvent);
		CL(Finish, cmdQueue);

		if (rank == 0) {
			cl_ulong timeStart, timeEnd;
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "Last Kerneltime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";

			uint64_t* R_ra = reflectancePerSimulation[simIndex];
			uint64_t* T_ra = transmissionPerSimulation[simIndex];
			uint64_t* A_rz = absorptionPerSimulation[simIndex];

			Write_Simulation_Results(A_rz, T_ra, R_ra, &simulations[simIndex], timeEnd - timeStart);

			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "Transfertime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";
		}

		CL(ReleaseEvent, kernelEvent); CL(ReleaseEvent, reflectanceTransferEvent);
	}
	freeResources();
}