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

using Weight = uint64_t;
using Buffer = void*;

// Representation of input file data
static SimulationStruct* simData = 0;
static int simCount = 0;

// 2D array with pointers to input buffers per simulation
static Buffer* input; // simCount x inBufCount
// 2D array with sizes of input buffers
static size_t* inBufSizes;
// Input buffer enumeration
enum InBufIndex { layerBuf, inBufCount };
// Get absolute index of input buffer
static int inBuf(int simIndex, InBufIndex inBuf) { return inBufCount * simCount + inBuf; }

// 2D array with pointers to output buffers per simulation
static Buffer* output; // simCount x outBufCount
// 2D array with sizes of ouput buffers
static size_t* outBufSizes;
// Output buffer enumeration
enum OutBufIndex { Rd, A, T, photonBuf, debugBuf, outBufCount };
// Get absolute index of output buffer
static int outBuf(int simIndex, OutBufIndex outBuf) { return outBufCount * simIndex + outBuf; }

// Triggers creation of debug buffer
static bool kernelDebugMode = false;

// Create new initial photon state
static PhotonState initPhotonState() {
	return {
		0.0f, 0.0f, 0.0f, // start position
		0.0f, 0.0f, 1.0f, // start direction
		1.0f, // start weight
		0 // start layer index
	};
}

// Fill layer buffer with data
static void initLayerBuf(int layerCount, Layer* layerBuf, SimulationStruct data) {
	for (int layerIndex = 1; layerIndex <= layerCount; layerIndex++) {
		layerBuf[layerIndex - 1] = {
			data.layers[layerIndex].mua,
			1.0f / data.layers[layerIndex].mutr - data.layers[layerIndex].mua,
			data.layers[layerIndex].g,
			data.layers[layerIndex].n,
			Boundary{data.layers[layerIndex].z_min, 0.0f, 0.0f, 1.0f},
			Boundary{data.layers[layerIndex].z_max, 0.0f, 0.0f, -1.0f},
		};
	}
}

// Interface for preparing kernel resources
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* mcmlOptions,
int* inputBufferCount, size_t* inputBufferSizes,
int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount) {
	// Read input file
	//TODO read with rank 0 and broadcast input data
	int ignoreA = std::string(kernelOptions).find("-D IGNORE_A") != std::string::npos ? 1 : 0;
	std::cout << "--- "<<mcmlOptions<<" --->" << std::endl;
	simCount = read_simulation_data(mcmlOptions, &simData, ignoreA);
	std::cout << "<--- "<<mcmlOptions<<" ---" << std::endl;

	// Ensure no bins can overflow
	for (int simIndex = 0; simIndex < simCount; simIndex++) {
		assert(simData[simIndex].number_of_photons <= 0xFFFFFFFFu);
	}

	// Alloc input buffers
	*inputBufferCount = simCount * inBufCount;
	assert(*inputBufferCount <= maxBufferCount);
	input = (Buffer*)malloc(simCount * inBufCount * sizeof(Buffer*));
	for (int simIndex = 0; simIndex < simCount; simIndex++) {
		// Layer buffer
		int layerCount = simData[simIndex].n_layers;
		size_t layerBufSize = layerCount * sizeof(Layer);
		input[inBuf(simIndex, layerBuf)] = malloc(layerBufSize);
		initLayerBuf(layerCount, (Layer*)input[inBuf(simIndex, layerBuf)], simData[simIndex]);
		inputBufferSizes[inBuf(simIndex, layerBuf)] = layerBufSize;
	}
	inBufSizes = inputBufferSizes;

	// Alloc output buffers
	kernelDebugMode = std::string(kernelOptions).find("-D DEBUG") != std::string::npos ? 1 : 0;
	*outputBufferCount = kernelDebugMode ? simCount * outBufCount : simCount * (outBufCount - 1);
	assert(*outputBufferCount <= maxBufferCount);
	output = (Buffer*)malloc(simCount * outBufCount * sizeof(Buffer*));
	for (int simIndex = 0; simIndex < simCount; simIndex++) {
		// Reflectance buffer
		int radialBinCount = simData[simIndex].det.nr;
		int angularBinCount = simData[simIndex].det.na;
		size_t reflectanceBufSize = radialBinCount * angularBinCount * sizeof(Weight);
		output[outBuf(simIndex, Rd)] = malloc(reflectanceBufSize);
		outputBufferSizes[outBuf(simIndex, Rd)] = reflectanceBufSize;
		// Absorption buffer
		int depthBinCount = simData[simIndex].det.nz;
		size_t absorptionBufSize = radialBinCount * depthBinCount * sizeof(Weight);
		output[outBuf(simIndex, A)] = malloc(absorptionBufSize);
		outputBufferSizes[outBuf(simIndex, A)] = absorptionBufSize;
		// Transmission buffer
		size_t transmissionBufSize = reflectanceBufSize;
		output[outBuf(simIndex, T)] = malloc(transmissionBufSize);
		outputBufferSizes[outBuf(simIndex, T)] = transmissionBufSize;
		// Photon state buffer
		size_t photonBufSize = totalThreadCount * sizeof(PhotonState);
		output[outBuf(simIndex, photonBuf)] = malloc(photonBufSize);
		outputBufferSizes[outBuf(simIndex, photonBuf)] = photonBufSize;
		// Debug buffer
		if (kernelDebugMode) {
			size_t debugBufSize = 2048;
			output[outBuf(simIndex, debugBuf)] = malloc(debugBufSize);
			outputBufferSizes[outBuf(simIndex, debugBuf)] = debugBufSize;
		}
	}
	outBufSizes = outputBufferSizes;
}

static void freeResources() {
	for (int i = 0; i < simCount; i++) {
		free(simData[i].layers);
		free(input[inBuf(i, layerBuf)]);
		free(output[outBuf(i, Rd)]);
		free(output[outBuf(i, A)]);
		free(output[outBuf(i, T)]);
		free(output[outBuf(i, photonBuf)]);
		if (kernelDebugMode) free(output[outBuf(i, debugBuf)]);
	}
	free(input);
	free(output);
	free(simData);
}

// Handle kernel debug output
static bool handleDebugOutput(int simIndex) {
	const char* error = "error";
	bool isError = false;
	int j = 0;
	char* debugBuffer = (char*)output[outBuf(simIndex, debugBuf)];
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

// Return name of the kernel function
const char* getCLKernelName() {
	return "mcml";
}

// Interface for running the kernel
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputGPU, cl_mem* outputGPU,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		// Upload layers
		CL(EnqueueWriteBuffer, cmdQueue, inputGPU[inBuf(simIndex, layerBuf)], CL_FALSE, 0,
			inBufSizes[inBuf(simIndex, layerBuf)], input[inBuf(simIndex, layerBuf)], 0, NULL, NULL);

		// First and last layer have only n initialized and represent outer media
		float nAbove = simData[simIndex].layers[0].n;
		float nBelow = simData[simIndex].layers[simData[simIndex].n_layers + 1].n;

		// Reflectance (specular):
		// percentage of light leaving at surface without any interaction
		// using Fesnel approximation by Schlick (no incident angle, no polarization)
		// Q: why are the ^2 different than in Schlicks approximation?
		float nDiff = nAbove - ((Layer*)input[inBuf(simIndex, layerBuf)])->n;
		float nSum = nAbove + ((Layer*)input[inBuf(simIndex, layerBuf)])->n;
		float R_specular = (nDiff * nDiff) / (nSum * nSum);
		
		// Init photon states
		for (int i = 0; i < totalThreadCount; i++) {
			PhotonState newState = initPhotonState();
			newState.weight -= R_specular;
			PhotonState* states = (PhotonState*)output[outBuf(simIndex, photonBuf)];
			states[i] = newState;
		}

		// Init accumulation buffers with zeros
		int radialBinCount = simData[simIndex].det.nr;
		int angularBinCount = simData[simIndex].det.na;
		for (int i = 0; i < radialBinCount * angularBinCount; i++) {
			Weight* weights = (Weight*)output[outBuf(simIndex, Rd)];
			weights[i] = 0;
		}

		// Upload reflectance buffer
		CL(EnqueueWriteBuffer, cmdQueue, outputGPU[outBuf(simIndex, Rd)], CL_FALSE, 0,
			outBufSizes[outBuf(simIndex, Rd)], output[outBuf(simIndex, Rd)], 0, NULL, NULL);

		{ // Set arguments
			float radialBinCentimeters = simData[simIndex].det.dr;
			int argCount = 0;
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &inputGPU[inBuf(simIndex, layerBuf)]);
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &simData[simIndex].n_layers);
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount);
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount);
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters);
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputGPU[outBuf(simIndex, Rd)]);
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputGPU[outBuf(simIndex, photonBuf)]);
			if (kernelDebugMode) {
				CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &outputGPU[outBuf(simIndex, debugBuf)]);
			}
		}

		// Run kernel with optimal thread count as long as targeted number of photons allows it
		//TODO compare perf against CUDAMCML, which does not wait for longest simulating thread,
		// but instead launches a new photon directly from a thread that would terminate,
		// causing additional tracking overhead.
		cl_event kernelEvent, reflectanceTransferEvent;
		uint32_t targetPhotonCount = simData[simIndex].number_of_photons;
		uint32_t finishedPhotonCount = 0;
		std::cout << std::endl;
		while (finishedPhotonCount < targetPhotonCount) { // stop when target reached
			// Upload photon states
			//TODO since buffer updates are sparse, map could be faster than write in whole
			CL(EnqueueWriteBuffer, cmdQueue, outputGPU[outBuf(simIndex, photonBuf)], CL_FALSE, 0,
				outBufSizes[outBuf(simIndex, photonBuf)], output[outBuf(simIndex, photonBuf)], 0, NULL, NULL);
			// Run a batch of photons
			size_t remainingPhotonCount = targetPhotonCount - finishedPhotonCount;
			if (remainingPhotonCount > totalThreadCount) {
				CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL,
					&totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
			} else {
				CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL,
					&remainingPhotonCount, NULL, 0, NULL, &kernelEvent);
			}
			// Download photon states
			CL(EnqueueReadBuffer, cmdQueue, outputGPU[outBuf(simIndex, photonBuf)], CL_FALSE, 0,
				outBufSizes[outBuf(simIndex, photonBuf)], output[outBuf(simIndex, photonBuf)], 0, NULL, NULL);
			if (kernelDebugMode) {
				CL(EnqueueReadBuffer, cmdQueue, outputGPU[outBuf(simIndex, debugBuf)], CL_FALSE, 0,
					2048, output[outBuf(simIndex, debugBuf)], 0, NULL, NULL);
			}
			// Wait for async commands to finish
			CL(Finish, cmdQueue);
			if (kernelDebugMode) {
				assert(!handleDebugOutput(simIndex));
			}
			// Check for dead photons
			for (int i = 0; i < totalThreadCount; i++) {
				PhotonState* states = (PhotonState*)output[outBuf(simIndex, photonBuf)];
				if (states[i].weight == 0) {
					finishedPhotonCount++;
					PhotonState newState = initPhotonState();
					newState.weight -= R_specular;
					states[i] = newState;
				}
			}
			if (rank == 0) {
				std::cout << '\r' << "Photons terminated: " << finishedPhotonCount << "/" << targetPhotonCount << std::flush;
			}
		}
		std::cout << std::endl;

		// Download reflectance
		CL(EnqueueReadBuffer, cmdQueue, outputGPU[outBuf(simIndex, Rd)], CL_FALSE, 0,
			outBufSizes[outBuf(simIndex, Rd)], output[outBuf(simIndex, Rd)], 0, NULL, &reflectanceTransferEvent);
		CL(Finish, cmdQueue);

		// Write output
		if (rank == 0) {
			cl_ulong timeStart, timeEnd;
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "Last Kerneltime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";

			Write_Simulation_Results(
				(Weight*)output[outBuf(simIndex, A)],
				(Weight*)output[outBuf(simIndex, T)],
				(Weight*)output[outBuf(simIndex, Rd)],
				&simData[simIndex],
				timeEnd - timeStart);

			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &timeStart, NULL);
			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &timeEnd, NULL);
			std::cout << "Transfertime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";
		}

		CL(ReleaseEvent, kernelEvent); CL(ReleaseEvent, reflectanceTransferEvent);
	}
	freeResources();
}
