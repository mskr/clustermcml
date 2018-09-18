#include <iostream> // std::cout
#include <string> // std::string
#include <stdint.h> // uint32_t, uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

#ifndef CL2CPP
#define DEBUG
#include "clcheck.h"
#endif

#include "BoundaryLayerPhoton.h"

#ifdef CL2CPP
#define cl_int int
#define cl_context int
#define cl_command_queue int
#define cl_kernel int
#define cl_mem int
// link or include?
// void mcml(float nAbove, float nBelow, struct Layer* layers, int layerCount,
// 	int size_r, int size_a, float delta_r, 
// 	volatile uint64_t* R_ra,
// 	struct PhotonState* photonStates);
#include "kernel.c"
#endif

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


void allocCLKernelResources(size_t totalThreadCount, const char* kernelOptions, char* mcmlOptions,
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

	assert(*inputBufferCount <= maxBufferCount);
	assert(*outputBufferCount <= maxBufferCount);

	for (int i = 0; i < simCount; i++) {

		assert(simulations[i].number_of_photons <= 0xFFFFFFFFu); // ensures no bins can overflow

		// Layers buffer
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

		// Reflectance buffer
		int radialBinCount = simulations[i].det.nr;
		int angularBinCount = simulations[i].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);
		uint64_t* R_ra = (uint64_t*)malloc(reflectanceBufferSize);
		reflectancePerSimulation[i] = R_ra;
		outputBufferSizes[i] = reflectanceBufferSize;

		// Transmission buffer
		size_t transmissionBufferSize = reflectanceBufferSize;
		uint64_t* T_ra = (uint64_t*)malloc(transmissionBufferSize);
		transmissionPerSimulation[i] = T_ra;

		// Absorption buffer
		int depthBinCount = simulations[i].det.nz;
		size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(uint64_t);
		uint64_t* A_rz = (uint64_t*)malloc(absorptionBufferSize);
		absorptionPerSimulation[i] = A_rz;

		// Photon states buffer
		photonStatesPerSimulation[i] = (PhotonState*)malloc(totalThreadCount * sizeof(PhotonState));
		outputBufferSizes[simCount + i] = totalThreadCount * sizeof(PhotonState);
	}

	// Debug buffer
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
		free(absorptionPerSimulation[i]);
		free(transmissionPerSimulation[i]);
		free(reflectancePerSimulation[i]);
		free(layersPerSimulation[i]);
		free(simulations[i].layers);
	}
	free(absorptionPerSimulation);
	free(transmissionPerSimulation);
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
		float nDiff = nAbove - simulations[simIndex].layers[1].n;
		float nSum = nAbove + simulations[simIndex].layers[1].n;
		float R_specular = (nDiff * nDiff) / (nSum * nSum);
		
		for (int i = 0; i < totalThreadCount; i++) {
			PhotonState newState = createNewPhotonState();
			newState.weight -= R_specular;
			photonStatesPerSimulation[simIndex][i] = newState;
		}

		int radialBinCount = simulations[simIndex].det.nr;
		float radialBinCentimeters = simulations[simIndex].det.dr;
		int angularBinCount = simulations[simIndex].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);

		// Init accumulation buffers with zeros
		for (int i = 0; i < radialBinCount * angularBinCount; i++) {
			reflectancePerSimulation[simIndex][i] = 0;
		}

		size_t photonStateBufferSize = totalThreadCount * sizeof(PhotonState);

		// Run kernel with optimal thread count as long as targeted number of photons allows it
		//TODO compare perf against CUDAMCML, which does not wait for longest simulating thread,
		// but instead launches a new photon directly from a thread that would terminate,
		// causing additional tracking overhead.
		uint32_t targetPhotonCount = simulations[simIndex].number_of_photons;
		uint32_t finishedPhotonCount = 0;
		std::cout << std::endl;
		while (finishedPhotonCount < targetPhotonCount) {

			mcml(nAbove, nBelow, layersPerSimulation[simIndex], simulations[simIndex].n_layers, // input
				radialBinCount, angularBinCount, radialBinCentimeters, reflectancePerSimulation[simIndex], // output
				photonStatesPerSimulation[simIndex]);// intermediate buffer

			// Check for dead photons
			for (int i = 0; i < totalThreadCount; i++) {
				if (photonStatesPerSimulation[simIndex][i].weight == 0) {
					finishedPhotonCount++;
					PhotonState newState = createNewPhotonState();
					newState.weight -= R_specular;
					photonStatesPerSimulation[simIndex][i] = newState;
				}
			}
			if (rank == 0) {
				std::cout << '\r' << "Photons terminated: " << finishedPhotonCount << "/" << targetPhotonCount << std::flush;
			}
		}
		std::cout << std::endl;

		// Write output
		if (rank == 0) {

			uint64_t* R_ra = reflectancePerSimulation[simIndex];
			uint64_t* T_ra = transmissionPerSimulation[simIndex];
			uint64_t* A_rz = absorptionPerSimulation[simIndex];

			Write_Simulation_Results(A_rz, T_ra, R_ra, &simulations[simIndex], 0);
		}
	}
	freeResources();
}

// The changes to run cl kernel as normal c++ function were only marginal
//TODO integrate this as option into main.cpp

int main(int nargs, char* args[]) {
	assert(nargs == 2);
	char* filename = args[1];
	std::cout << filename << std::endl;
	{
		int inputBufferCount = 0, outputBufferCount = 0;
		size_t inputBufferSizes[10], outputBufferSizes[10];
		cl_mem inputBuffers[10], outputBuffers[10];
		allocCLKernelResources(1, "", filename, &inputBufferCount, inputBufferSizes, &outputBufferCount, outputBufferSizes, 10);
	}

	runCLKernel((cl_context)0, (cl_command_queue)0, (cl_kernel)0, 0, 0, 1, 1, 1, 0);
	return 0;
}