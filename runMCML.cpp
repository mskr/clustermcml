#include <iostream> // std::cout
#include <string> // std::string
#include <stdint.h> // uint32_t, uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

#define DEBUG
#include "clcheck.h"

//TODO share random functions with kernel via header
#ifdef CL2CPU
// Need declarations when linking with kernel file as CPU code
uint32_t wang_hash(uint32_t seed);
void mcml(float nAbove, float nBelow, struct Layer* layers, int layerCount,
	int size_r, int size_a, int size_z, float delta_r,  float delta_z,
	volatile uint64_t* R_ra, volatile uint64_t* T_ra, volatile uint64_t* A_rz,
	struct PhotonState* photonStates);
#else
// Need definition when kernel is GPU code
uint32_t wang_hash(uint32_t seed) {
	seed = (seed ^ 61) ^ (seed >> 16);
	seed *= 9;
	seed = seed ^ (seed >> 4);
	seed *= 0x27d4eb2d;
	seed = seed ^ (seed >> 15);
	return seed;
}
#endif

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
	unsigned int rngState;
	unsigned int isDead;
};


static PhotonState createNewPhotonState() {
	return {
		0.0f, 0.0f, 0.0f, // start position
		0.0f, 0.0f, 1.0f, // start direction
		1.0f, // start weight
		0, // start layer index
		0, // dummy rng state
		0 // alive
	};
}


const char* getCLKernelName() {
	return "mcml";
}


static int simCount = 0;
static SimulationStruct* simulations = 0;
static Layer** layersPerSimulation = 0;
static uint64_t** reflectancePerSimulation = 0;
static uint64_t** transmissionPerSimulation = 0;
static uint64_t** absorptionPerSimulation = 0;

static PhotonState* stateBuffer = 0;
static char* debugBuffer = 0;

/**
* Alloc host buffers
* Load input data into host buffers
* Report size information
*/
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* mcmlOptions,
int* inputBufferCount, size_t* inputBufferSizes,
int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount) {

	//TODO read with rank 0 and broadcast input data
	int ignoreA = std::string(kernelOptions).find("-D IGNORE_A") != std::string::npos ? 1 : 0;
	std::cout << "--- "<<mcmlOptions<<" --->" << std::endl;
	simCount = read_simulation_data(mcmlOptions, &simulations, ignoreA);
	strcpy(simulations[0].outp_filename, "own.mco");
	std::cout << "<--- "<<mcmlOptions<<" ---" << std::endl;

	layersPerSimulation = (Layer**)malloc(simCount * sizeof(Layer*));
	reflectancePerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));
	transmissionPerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));
	absorptionPerSimulation = (uint64_t**)malloc(simCount * sizeof(uint64_t*));

	// Photon states buffer
	stateBuffer = (PhotonState*)malloc(totalThreadCount * sizeof(PhotonState));
	outputBufferSizes[0] = totalThreadCount * sizeof(PhotonState);

	*inputBufferCount = simCount;
	*outputBufferCount = 1 + 3 * simCount;

	assert(*inputBufferCount <= maxBufferCount);
	assert(*outputBufferCount <= maxBufferCount);

	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		assert(simulations[simIndex].number_of_photons <= 0xFFFFFFFFu); // ensures no bins can overflow

		// Layers buffer
		int layerCount = simulations[simIndex].n_layers;
		inputBufferSizes[simIndex] = layerCount * sizeof(Layer);
		Layer* layers = (Layer*)malloc(layerCount * sizeof(Layer));
		layersPerSimulation[simIndex] = layers;
		for (int j = 1; j <= layerCount; j++) {
			layers[j - 1] = {
				simulations[simIndex].layers[j].mua,
				1.0f / simulations[simIndex].layers[j].mutr - simulations[simIndex].layers[j].mua,
				simulations[simIndex].layers[j].g,
				simulations[simIndex].layers[j].n,
				Boundary{simulations[simIndex].layers[j].z_min, 0.0f, 0.0f, 1.0f},
				Boundary{simulations[simIndex].layers[j].z_max, 0.0f, 0.0f, -1.0f},
			};
		}

		// Reflectance buffer
		int radialBinCount = simulations[simIndex].det.nr;
		int angularBinCount = simulations[simIndex].det.na;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);
		uint64_t* R_ra = (uint64_t*)malloc(reflectanceBufferSize);
		reflectancePerSimulation[simIndex] = R_ra;
		outputBufferSizes[1 + simIndex] = reflectanceBufferSize;

		// Transmission buffer
		size_t transmissionBufferSize = reflectanceBufferSize;
		uint64_t* T_ra = (uint64_t*)malloc(transmissionBufferSize);
		transmissionPerSimulation[simIndex] = T_ra;
		outputBufferSizes[1 + simCount + simIndex] = transmissionBufferSize;

		// Absorption buffer
		int depthBinCount = simulations[simIndex].det.nz;
		size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(uint64_t);
		uint64_t* A_rz = (uint64_t*)malloc(absorptionBufferSize);
		absorptionPerSimulation[simIndex] = A_rz;
		outputBufferSizes[1 + 2 * simCount + simIndex] = absorptionBufferSize;
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
		free(stateBuffer);
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
	std::cout << std::endl;
	//TODO print everything as hex view until reaching some unique end symbol
	int n = 2048/4;
	float sum = 0.0f;
	for (int k = 0; k < n; k++) { // print as floats
		float f = ((float*)debugBuffer)[k];
		sum += f;
		std::cout << f << " ";
	}
	std::cout << "sum=" << sum << std::endl;
	std::cout << "average=" << (sum / n) << std::endl;
	std::cout << std::endl;
	return true;
}

/**
* Upload input data to device buffers
* Run kernel
* Download output data to host buffers
* Write output
*/
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		// Upload layers
		CL(EnqueueWriteBuffer, cmdQueue, inputBuffers[simIndex], CL_FALSE, 0,
			simulations[simIndex].n_layers * sizeof(Layer), layersPerSimulation[simIndex], 0, NULL, NULL);

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
		
		// Init photon states
		for (int i = 0; i < totalThreadCount; i++) {
			PhotonState newState = createNewPhotonState();
			newState.weight -= R_specular;
			stateBuffer[i] = newState;
		}

		// Get RAT buffer info
		float radialBinCentimeters = simulations[simIndex].det.dr;
		float depthBinCentimeters = simulations[simIndex].det.dz;
		int radialBinCount = simulations[simIndex].det.nr;
		int angularBinCount = simulations[simIndex].det.na;
		int depthBinCount = simulations[simIndex].det.nz;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(uint64_t);
		size_t transmissionBufferSize = reflectanceBufferSize;
		size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(uint64_t);

		// Init RAT buffers with zeros
		for (int i = 0; i < radialBinCount * angularBinCount; i++)
			reflectancePerSimulation[simIndex][i] = transmissionPerSimulation[simIndex][i] = 0;
		for (int i = 0; i < radialBinCount * depthBinCount; i++)
			absorptionPerSimulation[simIndex][i] = 0;

		// Upload RAT buffers
		CL(EnqueueWriteBuffer, cmdQueue, outputBuffers[1 + simIndex], CL_FALSE, 0,
			reflectanceBufferSize, reflectancePerSimulation[simIndex], 0, NULL, NULL);
		CL(EnqueueWriteBuffer, cmdQueue, outputBuffers[1 + simCount + simIndex], CL_FALSE, 0,
			transmissionBufferSize, transmissionPerSimulation[simIndex], 0, NULL, NULL);
		CL(EnqueueWriteBuffer, cmdQueue, outputBuffers[1 + 2 * simCount + simIndex], CL_FALSE, 0,
			absorptionBufferSize, absorptionPerSimulation[simIndex], 0, NULL, NULL);

		{ // Set arguments
			int argCount = 0;
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &inputBuffers[simIndex]); // layers
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &simulations[simIndex].n_layers); // layerCount
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount); // size_r
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount); // size_a
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &depthBinCount); // size_z
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters); // delta_r
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &depthBinCentimeters); // delta_z
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[1 + simIndex]); // R
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[1 + simCount + simIndex]); // T
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[1 + 2 * simCount + simIndex]); // A
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[0]); // Photon states
			if (debugBuffer) {
				CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputBuffers[1 + 3 * simCount]);
			}
		}

		cl_event kernelEvent, reflectanceTransferEvent;

		size_t photonStateBufferSize = totalThreadCount * sizeof(PhotonState);

		// Run kernel with optimal thread count as long as targeted number of photons allows it
		//TODO compare perf against CUDAMCML, which does not wait for longest simulating thread,
		// but instead launches a new photon directly from a thread that would terminate,
		// causing additional tracking overhead.
		uint32_t targetPhotonCount = simulations[simIndex].number_of_photons;
		uint32_t finishedPhotonCount = 0;
		std::cout << std::endl;
		while (finishedPhotonCount < targetPhotonCount) { // stop when target reached

			// Upload photon states
			//TODO since buffer updates are sparse, map could be faster than write in whole
			CL(EnqueueWriteBuffer, cmdQueue, outputBuffers[0], CL_FALSE, 0,
				photonStateBufferSize, stateBuffer, 0, NULL, NULL);

			// Run a batch of photons
			CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
			#ifdef CL2CPU
				mcml(nAbove, nBelow, layersPerSimulation[simIndex], simulations[simIndex].n_layers, // input
					radialBinCount, angularBinCount, depthBinCount, radialBinCentimeters, depthBinCentimeters
					reflectancePerSimulation[simIndex], transmissionPerSimulation[simIndex] absorptionPerSimulation[simIndex], // output
					stateBuffer);// intermediate buffer
			#endif

			// Download photon states
			CL(EnqueueReadBuffer, cmdQueue, outputBuffers[0], CL_FALSE, 0,
				photonStateBufferSize, stateBuffer, 0, NULL, NULL);
			if (debugBuffer) {
				CL(EnqueueReadBuffer, cmdQueue, outputBuffers[1 + 3 * simCount], CL_FALSE, 0, 2048, debugBuffer, 0, NULL, NULL);
			}

			// Wait for async commands to finish
			CL(Finish, cmdQueue);
			if (debugBuffer) {
				if (handleDebugOutput()) exit(1);
			}

			// Check for finished photons
			for (int i = 0; i < totalThreadCount; i++) {
				if (stateBuffer[i].weight == 0 && !(stateBuffer[i].isDead)) {
					finishedPhotonCount++;
					// Launch new only if next round cannot overachieve
					if (finishedPhotonCount + totalThreadCount <= targetPhotonCount) {
						PhotonState newState = createNewPhotonState();
						newState.weight -= R_specular;
						newState.rngState = wang_hash(finishedPhotonCount + totalThreadCount);
						stateBuffer[i] = newState;
					} else {
						stateBuffer[i].isDead = 1;
					}
				}
			}

			// Print status
			if (rank == 0) {
				std::cout << '\r' << "Photons terminated: " << finishedPhotonCount << "/" << targetPhotonCount << std::flush;
			}
		}
		std::cout << std::endl;

		// Download RAT
		CL(EnqueueReadBuffer, cmdQueue, outputBuffers[1 + simIndex], CL_FALSE, 0,
			reflectanceBufferSize, reflectancePerSimulation[simIndex], 0, NULL, &reflectanceTransferEvent);
		CL(EnqueueReadBuffer, cmdQueue, outputBuffers[1 + simCount + simIndex], CL_FALSE, 0,
			transmissionBufferSize, transmissionPerSimulation[simIndex], 0, NULL, NULL);
		CL(EnqueueReadBuffer, cmdQueue, outputBuffers[1 + 2 * simCount + simIndex], CL_FALSE, 0,
			absorptionBufferSize, absorptionPerSimulation[simIndex], 0, NULL, NULL);
		CL(Finish, cmdQueue);

		// Write output
		if (rank == 0) {
			cl_ulong timeStart = 0, timeEnd = 0;
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
