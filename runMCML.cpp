#include <string.h> // strstr, memcpy
#include <stdint.h> // uint32_t, uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

#define DEBUG
#include "clcheck.h" // CL macro
#include "mpicheck.h" // MPI macro

#include "Log.h" // out stream
#include "randomlib.h" // wang_hash, rand_xorshift

#include "clmem.h" // CLMALLOC_INPUT, CLMALLOC_OUTPUT, getCLHostPointer, CLMEM_ACCESS_ARRAY, CLMEM_ACCESS_AOS

// Data structures
#include "Boundary.h"
#include "Layer.h"
#include "PhotonTracker.h"

#include "CUDAMCML.h"
#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

typedef uint64_t Weight;
#define MPI_WEIGHT_T MPI_UINT64_T

// Declaration of kernel function to be able to use it from CPU code
void mcml(float nAbove, float nBelow, struct Layer* layers, int layerCount,
int size_r, int size_a, int size_z, float delta_r,  float delta_z,
volatile Weight* R_ra, volatile Weight* T_ra, volatile Weight* A_rz,
struct PhotonTracker* photonStates);


static PhotonTracker createNewPhotonTracker() {
	return {
		0.0f, 0.0f, 0.0f, // start position
		0.0f, 0.0f, 1.0f, // start direction
		1.0f, // start weight
		0, // start layer index
		0, // dummy rng state
		0 // alive
	};
}


static MPI_Datatype createMPISimulationStruct() {
	// Describing the following data layout:
	// uint32_t number_of_photons;
	// uint32_t ignoreAdetection;
	// uint32_t n_layers;
	// uint32_t start_weight;
	// uint32_t begin,end; 		// mci file position offsets
	// char outp_filename[STR_LEN];
	// char inp_filename[STR_LEN];
	// char AorB, padding[7]; 	// enforce 8 byte alignment
	// struct DetStruct {
	// 	float dr;				// Detection grid resolution, r-direction [cm]
	// 	float dz;				// Detection grid resolution, z-direction [cm]
	// 	uint32_t na;			// Number of grid elements in angular-direction [-]
	// 	uint32_t nr;			// Number of grid elements in r-direction
	// 	uint32_t nz;			// Number of grid elements in z-direction
	// 	uint32_t padding; 		// enforce 8 byte alignment
	// } det;
	// LayerStruct* layers;
	// Boundary* boundaries;
	int blockLengths[5] = {6,STR_LEN*2+8,2,4,2*sizeof(void*)};
	int offsets[5]; int sum = 0; offsets[0] = 0;
	offsets[1] = (sum += sizeof(uint32_t) * 6);
	offsets[2] = (sum += STR_LEN * 2 + 8);
	offsets[3] = (sum += sizeof(float) * 2);
	offsets[4] = (sum += sizeof(uint32_t) * 4);
	MPI_Datatype types[5] = {MPI_UINT32_T,MPI_CHAR,MPI_FLOAT,MPI_UINT32_T,MPI_CHAR};
	MPI_Datatype simStruct;
	MPI(Type_create_struct, 5, blockLengths, offsets, types, &simStruct);
	MPI(Type_commit, &simStruct);
	return simStruct;
}


static MPI_Datatype createMPILayerStruct() {
	int blockLengths[6] = {1,1,1,1,1,1};
	int offsets[6]; for (int k = 0; k < 6; k++) offsets[k] = k * sizeof(float);
	MPI_Datatype types[6] = {MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT};
	MPI_Datatype layerStruct;
	MPI(Type_create_struct, 6, blockLengths, offsets, types, &layerStruct);
	MPI(Type_commit, &layerStruct);
	return layerStruct;
}


static int simCount = 0;
static SimulationStruct* simulations = 0;
static cl_mem* layersPerSimulation = 0;
static cl_mem* boundariesPerSimulation = 0;
static cl_mem* reflectancePerSimulation = 0;
static cl_mem* transmissionPerSimulation = 0;
static cl_mem* absorptionPerSimulation = 0;

static cl_mem stateBuffer = 0;
static cl_mem debugBuffer = 0;

/**
* Alloc host buffers
* Load input data into host buffers
* Report size information
*/
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* mcmlOptions, int rank) {

	// Read input file with process 0
	if (rank == 0) {
		int ignoreA = strstr(kernelOptions, "-D IGNORE_A") != NULL ? 1 : 0;
		out << "--- "<<mcmlOptions<<" --->" << '\n';
		simCount = read_simulation_data(mcmlOptions, &simulations, ignoreA, 1);
		out << "<--- "<<mcmlOptions<<" ---" << '\n';
	}

	// Broadcast the complex simulation struct piece by piece
	// (could optimize communication overhead with MPI pack/unpack:
	// https://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node84.htm
	// https://stackoverflow.com/a/32487093)
	MPI(Bcast, &simCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank != 0)
		simulations = (SimulationStruct*)malloc(simCount * sizeof(SimulationStruct));
	MPI_Datatype mpiSimStruct = createMPISimulationStruct();
	int mpiSimSize = 0;
	MPI(Type_size, mpiSimStruct, &mpiSimSize);
	assert(mpiSimSize == sizeof(SimulationStruct));
	MPI(Bcast, simulations, simCount, mpiSimStruct, 0, MPI_COMM_WORLD);
	MPI_Datatype mpiLayerStruct = createMPILayerStruct();
	int mpiLayerSize = 0;
	MPI(Type_size, mpiLayerStruct, &mpiLayerSize);
	assert(mpiLayerSize == sizeof(LayerStruct));
	for (int i = 0; i < simCount; i++) {
		if (rank != 0)
			simulations[i].layers = (LayerStruct*)malloc((simulations[i].n_layers + 2) * sizeof(LayerStruct));
		MPI(Bcast, simulations[i].layers, simulations[i].n_layers + 2, mpiLayerStruct, 0, MPI_COMM_WORLD);
	}

	// Allocate per-simulation arrays
	layersPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	reflectancePerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	transmissionPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	absorptionPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));

	// Photon states buffer
	stateBuffer = CLMALLOC_OUTPUT(totalThreadCount, PhotonTracker);

	// Allocate memory for each simulation
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		// ensure no bins can overflow
		assert(simulations[simIndex].number_of_photons <= 0xFFFFFFFFu);

		// Layers and boundaries buffers
		int layerCount = simulations[simIndex].n_layers;
		layersPerSimulation[simIndex] = CLMALLOC_INPUT(layerCount, Layer);
		boundariesPerSimulation[simIndex] = CLMALLOC_INPUT(layerCount+1, Boundary);
		for (int j = 1; j <= layerCount; j++) {
			// Boundary
			CLMEM_ACCESS_AOS(boundariesPerSimulation[simIndex], Boundary, j-1, z) = simulations[simIndex].layers[j].z_min;
			memcpy(CLMEM_ACCESS_AOS(boundariesPerSimulation[simIndex], Boundary, j-1, heightfield),
				simulations[simIndex].boundaries[j-1].heightfield, BOUNDARY_SAMPLES);
			// Layer
			CLMEM_ACCESS_AOS(layersPerSimulation[simIndex], Layer, j-1, absorbCoeff) = simulations[simIndex].layers[j].mua;
			CLMEM_ACCESS_AOS(layersPerSimulation[simIndex], Layer, j-1, scatterCoeff) = 
				1.0f / simulations[simIndex].layers[j].mutr - simulations[simIndex].layers[j].mua;
			CLMEM_ACCESS_AOS(layersPerSimulation[simIndex], Layer, j-1, g) = simulations[simIndex].layers[j].g;
			CLMEM_ACCESS_AOS(layersPerSimulation[simIndex], Layer, j-1, n) = simulations[simIndex].layers[j].n;
		}
		// One more boundary
		CLMEM_ACCESS_AOS(boundariesPerSimulation[simIndex], Boundary, layerCount, z) = simulations[simIndex].layers[layerCount-1].z_max;
		memcpy(CLMEM_ACCESS_AOS(boundariesPerSimulation[simIndex], Boundary, layerCount, heightfield), 
			simulations[simIndex].boundaries[layerCount].heightfield, BOUNDARY_SAMPLES);

		// Reflectance buffer
		int radialBinCount = simulations[simIndex].det.nr;
		int angularBinCount = simulations[simIndex].det.na;
		reflectancePerSimulation[simIndex] = CLMALLOC_OUTPUT(radialBinCount*angularBinCount, Weight);

		// Transmission buffer
		transmissionPerSimulation[simIndex] = CLMALLOC_OUTPUT(radialBinCount*angularBinCount, Weight);

		// Absorption buffer
		int depthBinCount = simulations[simIndex].det.nz;
		absorptionPerSimulation[simIndex] = CLMALLOC_OUTPUT(radialBinCount*depthBinCount, Weight);
	}

	// Debug buffer
	int debugMode = strstr(kernelOptions, "-D DEBUG") != NULL ? 1 : 0;
	if (debugMode) debugBuffer = CLMALLOC_OUTPUT(2048, char);
}


static void freeResources() {
	for (int i = 0; i < simCount; i++) {
		free(simulations[i].layers);
		free(simulations[i].boundaries);
	}
	free(absorptionPerSimulation);
	free(transmissionPerSimulation);
	free(reflectancePerSimulation);
	free(layersPerSimulation);
	free(simulations);
}


static bool handleDebugOutput() {
	out << '\n';
	//TODO print everything as hex view until reaching some unique end symbol
	int n = 2048/4;
	float sum = 0.0f;
	for (int k = 0; k < n; k++) { // print as floats
		float f = ((float*)getCLHostPointer(debugBuffer))[k];
		sum += f;
		out << f << " ";
	}
	out << "sum=" << sum << '\n';
	out << "average=" << (sum / n) << '\n';
	out << '\n';
	return true;
}


//TODO refactor into smaller functions

static void restartPhotons(PhotonTracker* stateBuffer, uint32_t finishCount, uint32_t targetCount) {

}


static void doOneGPURun(PhotonTracker* stateBuffer, uint32_t finishCount, uint32_t targetCount) {

}


static void doOneSimulation(float nAbove, float nBelow, Layer* layersHost, cl_mem layersDevice, int layerCount, size_t layerBufferSize) {

}

//TODO can we use async buffer transfers?


/**
* Upload input data to device buffers
* Run kernel
* Download output data to host buffers
* Write output
*/
void runCLKernel(cl_command_queue cmdQueue, cl_kernel kernel, size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		// Upload layers
		CL(EnqueueWriteBuffer, cmdQueue, layersPerSimulation[simIndex], CL_FALSE, 0,
			simulations[simIndex].n_layers * sizeof(Layer), getCLHostPointer(layersPerSimulation[simIndex]), 0, NULL, NULL);

		// first and last layer have only n initialized and represent outer media
		float nAbove = simulations[simIndex].layers[0].n;
		float nBelow = simulations[simIndex].layers[simulations[simIndex].n_layers + 1].n;

		// Reflectance (specular):
		// percentage of light leaving at surface without any interaction
		// using Fesnel approximation by Schlick (no incident angle, no polarization)
		// Q: why are the ^2 different than in Schlicks approximation?
		float nDiff = nAbove - CLMEM_ACCESS_AOS(layersPerSimulation[simIndex], Layer, 0, n);
		float nSum = nAbove + CLMEM_ACCESS_AOS(layersPerSimulation[simIndex], Layer, 0, n);
		float R_specular = (nDiff * nDiff) / (nSum * nSum);

		// Get RAT buffer info
		float radialBinCentimeters = simulations[simIndex].det.dr;
		float depthBinCentimeters = simulations[simIndex].det.dz;
		int radialBinCount = simulations[simIndex].det.nr;
		int angularBinCount = simulations[simIndex].det.na;
		int depthBinCount = simulations[simIndex].det.nz;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(Weight);
		size_t transmissionBufferSize = reflectanceBufferSize;
		size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(Weight);

		// Init RAT buffers with zeros
		for (int i = 0; i < radialBinCount * angularBinCount; i++) {
			CLMEM_ACCESS_ARRAY(reflectancePerSimulation[simIndex], Weight, i) = 0;
			CLMEM_ACCESS_ARRAY(transmissionPerSimulation[simIndex], Weight, i) = 0;
		}
		for (int i = 0; i < radialBinCount * depthBinCount; i++)
			CLMEM_ACCESS_ARRAY(absorptionPerSimulation[simIndex], Weight, i) = 0;

		// Upload RAT buffers
		CL(EnqueueWriteBuffer, cmdQueue, reflectancePerSimulation[simIndex], CL_FALSE, 0,
			reflectanceBufferSize, getCLHostPointer(reflectancePerSimulation[simIndex]), 0, NULL, NULL);
		CL(EnqueueWriteBuffer, cmdQueue, transmissionPerSimulation[simIndex], CL_FALSE, 0,
			transmissionBufferSize, getCLHostPointer(transmissionPerSimulation[simIndex]), 0, NULL, NULL);
		CL(EnqueueWriteBuffer, cmdQueue, absorptionPerSimulation[simIndex], CL_FALSE, 0,
			absorptionBufferSize, getCLHostPointer(absorptionPerSimulation[simIndex]), 0, NULL, NULL);

		{ // Set arguments
			int argCount = 0;
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &layersPerSimulation[simIndex]); // layers
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &simulations[simIndex].n_layers); // layerCount
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount); // size_r
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount); // size_a
			CL(SetKernelArg, kernel, argCount++, sizeof(int), &depthBinCount); // size_z
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters); // delta_r
			CL(SetKernelArg, kernel, argCount++, sizeof(float), &depthBinCentimeters); // delta_z
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &reflectancePerSimulation[simIndex]); // R
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &transmissionPerSimulation[simIndex]); // T
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &absorptionPerSimulation[simIndex]); // A
			CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &stateBuffer); // Photon states
			if (debugBuffer) CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &debugBuffer);
		}

		cl_event kernelEvent, reflectanceTransferEvent;

		size_t photonStateBufferSize = totalThreadCount * sizeof(PhotonTracker);

		uint32_t targetPhotonCount = simulations[simIndex].number_of_photons;

		// Distribute photons on processes
		uint32_t processPhotonCount = targetPhotonCount / processCount;
		if (rank == 0) { // rank 0 takes the remainder on top
			processPhotonCount += (targetPhotonCount % processCount);
		}
		
		// Init photon states
		for (int i = 0; i < totalThreadCount; i++) {
			if (i < targetPhotonCount) {
				PhotonTracker newState = createNewPhotonTracker();
				newState.weight -= R_specular;
				CLMEM_ACCESS_ARRAY(stateBuffer, PhotonTracker, i) = newState;
			} else {
				CLMEM_ACCESS_AOS(stateBuffer, PhotonTracker, i, isDead) = 1;
			}
		}

		// Run kernel with optimal thread count as long as targeted number of photons allows it
		//TODO compare perf against CUDAMCML, which does not wait for longest simulating thread,
		// but instead launches a new photon directly from a thread that would terminate,
		// causing additional tracking overhead.
		uint32_t finishedPhotonCount = 0;
		if (rank == 0) out << '\n';
		while (finishedPhotonCount < processPhotonCount) {

			// Upload photon states
			//TODO since buffer updates are sparse, map could be faster than write in whole
			CL(EnqueueWriteBuffer, cmdQueue, stateBuffer, CL_FALSE, 0,
				photonStateBufferSize, getCLHostPointer(stateBuffer), 0, NULL, NULL);

			// Run a batch of photons
			CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
			#ifdef CL2CPU
				mcml(nAbove, nBelow, getCLHostPointer(layersPerSimulation[simIndex]), simulations[simIndex].n_layers, // input
					radialBinCount, angularBinCount, depthBinCount, radialBinCentimeters, depthBinCentimeters,
					getCLHostPointer(reflectancePerSimulation[simIndex]),
					getCLHostPointer(transmissionPerSimulation[simIndex]),
					getCLHostPointer(absorptionPerSimulation[simIndex]), // output
					getCLHostPointer(stateBuffer));// intermediate buffer
			#endif

			// Download photon states
			CL(EnqueueReadBuffer, cmdQueue, stateBuffer, CL_FALSE, 0,
				photonStateBufferSize, getCLHostPointer(stateBuffer), 0, NULL, NULL);
			if (debugBuffer) { // download debug buffer if in debug mode
				CL(EnqueueReadBuffer, cmdQueue, debugBuffer, CL_FALSE, 0, 2048, getCLHostPointer(debugBuffer), 0, NULL, NULL);
			}

			// Wait for async commands to finish
			CL(Finish, cmdQueue);
			if (debugBuffer) {
				if (handleDebugOutput()) exit(1);
			}

			// Check for finished photons
			for (int i = 0; i < totalThreadCount; i++) {
				if (CLMEM_ACCESS_AOS(stateBuffer, PhotonTracker, i, weight) == 0
				&& !CLMEM_ACCESS_AOS(stateBuffer, PhotonTracker, i, isDead)) {
					finishedPhotonCount++;
					// Launch new only if next round cannot overachieve
					if (finishedPhotonCount + totalThreadCount <= processPhotonCount) {
						PhotonTracker newState = createNewPhotonTracker();
						newState.weight -= R_specular;
						newState.rngState = wang_hash(finishedPhotonCount + totalThreadCount);
						CLMEM_ACCESS_ARRAY(stateBuffer, PhotonTracker, i) = newState;
					} else {
						CLMEM_ACCESS_AOS(stateBuffer, PhotonTracker, i, isDead) = 1;
					}
				}
			}

			uint32_t totalFinishedPhotonCount = 0;
			MPI(Reduce, &finishedPhotonCount, &totalFinishedPhotonCount, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);

			// Print status
			if (rank == 0) out << '\r' << "Photons terminated: " << totalFinishedPhotonCount << "/" << targetPhotonCount << Log::flush;
		}
		if (rank == 0) out << '\n';

		// Download RAT
		CL(EnqueueReadBuffer, cmdQueue, reflectancePerSimulation[simIndex], CL_FALSE, 0,
			reflectanceBufferSize, getCLHostPointer(reflectancePerSimulation[simIndex]), 0, NULL, &reflectanceTransferEvent);
		CL(EnqueueReadBuffer, cmdQueue, transmissionPerSimulation[simIndex], CL_FALSE, 0,
			transmissionBufferSize, getCLHostPointer(transmissionPerSimulation[simIndex]), 0, NULL, NULL);
		CL(EnqueueReadBuffer, cmdQueue, absorptionPerSimulation[simIndex], CL_FALSE, 0,
			absorptionBufferSize, getCLHostPointer(absorptionPerSimulation[simIndex]), 0, NULL, NULL);
		CL(Finish, cmdQueue);

		// Sum RAT buffers from all processes
		Weight* totalReflectance = (Weight*)malloc(reflectanceBufferSize);
		MPI(Reduce, getCLHostPointer(reflectancePerSimulation[simIndex]), totalReflectance, radialBinCount*angularBinCount, 
			MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
		Weight* totalTransmittance = (Weight*)malloc(transmissionBufferSize);
		MPI(Reduce, getCLHostPointer(transmissionPerSimulation[simIndex]), totalTransmittance, radialBinCount*angularBinCount,
			MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
		Weight* totalAbsorption = (Weight*)malloc(absorptionBufferSize);
		MPI(Reduce, getCLHostPointer(absorptionPerSimulation[simIndex]), totalAbsorption, radialBinCount*depthBinCount, 
			MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);

		// Write output
		if (rank == 0) {
			uint64_t timeStart = 0, timeEnd = 0;
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(uint64_t), &timeStart, NULL);
			CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(uint64_t), &timeEnd, NULL);
			out << "Last Kerneltime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";

			Weight* R_ra = totalReflectance;
			Weight* T_ra = totalTransmittance;
			Weight* A_rz = totalAbsorption;
			Write_Simulation_Results(A_rz, T_ra, R_ra, &simulations[simIndex], timeEnd - timeStart);

			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(uint64_t), &timeStart, NULL);
			CL(GetEventProfilingInfo, reflectanceTransferEvent, CL_PROFILING_COMMAND_END, sizeof(uint64_t), &timeEnd, NULL);
			out << "Transfertime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";
		}

		free(totalReflectance);
		free(totalTransmittance);
		free(totalAbsorption);

		CL(ReleaseEvent, kernelEvent); CL(ReleaseEvent, reflectanceTransferEvent);
	}
	freeResources();
}


const char* getCLKernelName() {
	return "mcml";
}