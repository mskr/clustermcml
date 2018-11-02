/*********************************************************************************
*
* This code sets up MCML and distributes work on processes and threads.
*
* It runs one or more simulations from 1 mci file and produces mco files for each.
*
*********************************************************************************/

#include <string.h> // strstr, memcpy
#include <stdint.h> // uint32_t, uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

#include "clmem.h" // CLMALLOC_INPUT, CLMALLOC_OUTPUT, CLMEM, CLMEM_ACCESS_ARRAY, CLMEM_ACCESS_AOS

// Data structures
#include "Boundary.h"
#include "Layer.h"
#include "PhotonTracker.h"

#include "randomlib.h" // wang_hash, rand_xorshift
#include "Log.h" // out stream
#define DEBUG
#include "clcheck.h" // CL macro
#include "mpicheck.h" // MPI macro

#include "CUDAMCML.h"
#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

typedef uint64_t Weight;
#define MPI_WEIGHT_T MPI_UINT64_T


static int simCount = 0;
static SimulationStruct* simulations = 0;
static cl_mem* layersPerSimulation = 0;
static cl_mem* boundariesPerSimulation = 0;
static cl_mem* reflectancePerSimulation = 0;
static cl_mem* transmissionPerSimulation = 0;
static cl_mem* absorptionPerSimulation = 0;
static cl_mem stateBuffer = 0;
static cl_mem debugBuffer = 0;


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


static void readMCIFile(char* name, bool ignoreA, bool explicitBoundaries, int* outSimCount) {
	out << "Following info was read from input file \"" << name << "\":\n";
	*outSimCount = read_simulation_data(name, &simulations, ignoreA?1:0, explicitBoundaries?1:0);
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
	// Describing the following data layout:
	// float z_min;		// Layer z_min [cm]
	// float z_max;		// Layer z_max [cm]
	// float mutr;			// Reciprocal mu_total [cm]
	// float mua;			// Absorption coefficient [1/cm]
	// float g;			// Anisotropy factor [-]
	// float n;			// Refractive index [-]
	int blockLengths[1] = {6};
	int offsets[1] = {0};
	MPI_Datatype types[1] = {MPI_FLOAT};
	MPI_Datatype layerStruct;
	MPI(Type_create_struct, 1, blockLengths, offsets, types, &layerStruct);
	MPI(Type_commit, &layerStruct);
	return layerStruct;
}


static MPI_Datatype createMPIBoundaryStruct() {
	int blockLengths[1] = {1 + BOUNDARY_SAMPLES};
	int offsets[1] = {0};
	MPI_Datatype types[1] = {MPI_FLOAT};
	MPI_Datatype boundaryStruct;
	MPI(Type_create_struct, 1, blockLengths, offsets, types, &boundaryStruct);
	MPI(Type_commit, &boundaryStruct);
	return boundaryStruct;
}


static void broadcastInputData(int rank) {
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

	MPI_Datatype mpiBoundaryStruct = createMPIBoundaryStruct();
	int mpiBoundarySize = 0;
	MPI(Type_size, mpiBoundaryStruct, &mpiBoundarySize);
	assert(mpiBoundarySize == sizeof(Boundary));
	for (int i = 0; i < simCount; i++) {
		if (rank != 0)
			simulations[i].boundaries = (Boundary*)malloc((simulations[i].n_layers + 1) * sizeof(Boundary));
		MPI(Bcast, simulations[i].boundaries, simulations[i].n_layers + 1, mpiBoundaryStruct, 0, MPI_COMM_WORLD);
	}
}


static void setupInputArrays(SimulationStruct sim, int simIndex) {
	// Alloc space for layers and boundaries
	int layerCount = sim.n_layers;
	layersPerSimulation[simIndex] = CLMALLOC_INPUT(layerCount, Layer);
	boundariesPerSimulation[simIndex] = CLMALLOC_INPUT(layerCount+1, Boundary);

	// Data structures as read from file are slightly restructured to be ready for GPU consumption
	for (int j = 1; j <= layerCount; j++) {
		// Boundary
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, z) = sim.layers[j].z_min;
		memcpy(CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield),
			sim.boundaries[j-1].heightfield, BOUNDARY_SAMPLES);
		// Layer
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, absorbCoeff) = sim.layers[j].mua;
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, scatterCoeff) = 
			1.0f / sim.layers[j].mutr - sim.layers[j].mua;
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, g) = sim.layers[j].g;
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, n) = sim.layers[j].n;
	}
	// One more boundary
	CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, z) = sim.layers[layerCount].z_max;
	memcpy(CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield), 
		sim.boundaries[layerCount].heightfield, BOUNDARY_SAMPLES);
}


static void allocOutputArrays(SimulationStruct sim, int simIndex) {
	// Reflectance buffer
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	reflectancePerSimulation[simIndex] = CLMALLOC_OUTPUT(radialBinCount*angularBinCount, Weight);
	// Transmission buffer
	transmissionPerSimulation[simIndex] = CLMALLOC_OUTPUT(radialBinCount*angularBinCount, Weight);
	// Absorption buffer
	int depthBinCount = sim.det.nz;
	absorptionPerSimulation[simIndex] = CLMALLOC_OUTPUT(radialBinCount*depthBinCount, Weight);
}


static void freeResources() {
	for (int i = 0; i < simCount; i++) {
		free(simulations[i].layers);
		free(simulations[i].boundaries);
	}
	free(absorptionPerSimulation);
	free(transmissionPerSimulation);
	free(reflectancePerSimulation);
	free(boundariesPerSimulation);
	free(layersPerSimulation);
	free(simulations);
}


static bool handleDebugOutput() {
	out << '\n';
	//TODO print everything as hex view until reaching some unique end symbol
	int n = 2048/4;
	float sum = 0.0f;
	for (int k = 0; k < n; k++) { // print as floats
		float f = ((float*)CLMEM(debugBuffer))[k];
		sum += f;
		out << f << " ";
	}
	out << "sum=" << sum << '\n';
	out << "average=" << (sum / n) << '\n';
	out << '\n';
	return true;
}


static void uploadInputArrays(cl_command_queue cmd, SimulationStruct sim, cl_mem layers, cl_mem boundaries) {
	// Do a blocking write to be safe (but maybe slow)
	//TODO can we use async buffer transfers?
	CL(EnqueueWriteBuffer, cmd, layers, CL_TRUE, 0, sim.n_layers*sizeof(Layer), CLMEM(layers), 0, NULL, NULL);
	CL(EnqueueWriteBuffer, cmd, boundaries, CL_TRUE, 0, sim.n_layers*sizeof(Layer), CLMEM(boundaries), 0, NULL, NULL);
}


static void initAndUploadOutputArrays(cl_command_queue cmd, SimulationStruct sim, cl_mem R, cl_mem A, cl_mem T) {
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	int depthBinCount = sim.det.nz;
	size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(Weight);
	size_t transmissionBufferSize = reflectanceBufferSize;
	size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(Weight);

	// Init with zeros
	for (int i = 0; i < radialBinCount * angularBinCount; i++) {
		CLMEM_ACCESS_ARRAY(CLMEM(R), Weight, i) = 0;
		CLMEM_ACCESS_ARRAY(CLMEM(T), Weight, i) = 0;
	}
	for (int i = 0; i < radialBinCount * depthBinCount; i++)
		CLMEM_ACCESS_ARRAY(CLMEM(A), Weight, i) = 0;
	CL(EnqueueWriteBuffer, cmd, R, CL_TRUE, 0, reflectanceBufferSize, CLMEM(R), 0, NULL, NULL);
	CL(EnqueueWriteBuffer, cmd, T, CL_TRUE, 0, transmissionBufferSize, CLMEM(T), 0, NULL, NULL);
	CL(EnqueueWriteBuffer, cmd, A, CL_TRUE, 0, absorptionBufferSize, CLMEM(A), 0, NULL, NULL);
}


static float computeFresnelReflectance(float n0, float n1) {
	// Reflectance (specular):
	// percentage of light leaving at surface without any interaction
	// using Fesnel approximation by Schlick (no incident angle, no polarization)
	// Q: why are the ^2 different than in Schlicks approximation?
	float nDiff = n0 - n1;
	float nSum = n0 + n1;
	return (nDiff * nDiff) / (nSum * nSum);
}


static void initPhotonStates(size_t totalThreadCount, float R_specular) {
	for (int i = 0; i < totalThreadCount; i++) {
		PhotonTracker newState = createNewPhotonTracker();
		newState.weight -= R_specular;
		CLMEM_ACCESS_ARRAY(CLMEM(stateBuffer), PhotonTracker, i) = newState;
	}
}


static void restartFinishedPhotons(uint32_t processPhotonCount, size_t totalThreadCount, uint32_t* outFinishCount, float R_specular) {
	// Check for finished photons
	for (int i = 0; i < totalThreadCount; i++) {
		if (CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) == 0
		&& !CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead)) {
			(*outFinishCount)++;
			// Launch new only if next round cannot overachieve
			if ((*outFinishCount) + totalThreadCount <= processPhotonCount) {
				PhotonTracker newState = createNewPhotonTracker();
				newState.weight -= R_specular;
				newState.rngState = wang_hash((*outFinishCount) + totalThreadCount);
				CLMEM_ACCESS_ARRAY(CLMEM(stateBuffer), PhotonTracker, i) = newState;
			} else {
				CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead) = 1;
			}
		}
	}
}


static void setKernelArguments(cl_kernel kernel, SimulationStruct sim, cl_mem layers, cl_mem boundaries, cl_mem R, cl_mem A, cl_mem T) {
	float nAbove = sim.layers[0].n; // first and last layer have only n initialized...
	float nBelow = sim.layers[sim.n_layers + 1].n; // ...and represent outer media
	float radialBinCentimeters = sim.det.dr;
	float depthBinCentimeters = sim.det.dz;
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	int depthBinCount = sim.det.nz;
	int argCount = 0;
	CL(SetKernelArg, kernel, argCount++, sizeof(float), &nAbove);
	CL(SetKernelArg, kernel, argCount++, sizeof(float), &nBelow);
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &layers); // layers
	CL(SetKernelArg, kernel, argCount++, sizeof(int), &sim.n_layers); // layerCount
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &boundaries); // boundaries
	CL(SetKernelArg, kernel, argCount++, sizeof(int), &radialBinCount); // size_r
	CL(SetKernelArg, kernel, argCount++, sizeof(int), &angularBinCount); // size_a
	CL(SetKernelArg, kernel, argCount++, sizeof(int), &depthBinCount); // size_z
	CL(SetKernelArg, kernel, argCount++, sizeof(float), &radialBinCentimeters); // delta_r
	CL(SetKernelArg, kernel, argCount++, sizeof(float), &depthBinCentimeters); // delta_z
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &R); // R
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &T); // T
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &A); // A
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &stateBuffer); // Photon states
	if (debugBuffer) CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &debugBuffer);
}


static void simulate(cl_command_queue cmd, cl_kernel kernel, size_t totalThreadCount, size_t simdThreadCount, 
uint32_t photonCount, uint32_t targetPhotonCount, int rank, float R_specular, cl_event kernelEvent) {

	size_t photonStateBufferSize = totalThreadCount * sizeof(PhotonTracker);

	// Run kernel with optimal thread count as long as targeted number of photons allows it
	//TODO compare perf against CUDAMCML, which does not wait for longest simulating thread,
	// but instead launches a new photon directly from a thread that would terminate,
	// causing additional tracking overhead.
	uint32_t finishedPhotonCount = 0;
	if (rank == 0) out << '\n';
	while (finishedPhotonCount < photonCount) {

		// Upload photon states
		//TODO since buffer updates are sparse, map could be faster than write in whole
		CL(EnqueueWriteBuffer, cmd, stateBuffer, CL_FALSE, 0,
			photonStateBufferSize, CLMEM(stateBuffer), 0, NULL, NULL);

		// Run a batch of photons
		CL(EnqueueNDRangeKernel, cmd, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);

		// Download photon states
		CL(EnqueueReadBuffer, cmd, stateBuffer, CL_FALSE, 0,
			photonStateBufferSize, CLMEM(stateBuffer), 0, NULL, NULL);
		if (debugBuffer) { // download debug buffer if in debug mode
			CL(EnqueueReadBuffer, cmd, debugBuffer, CL_FALSE, 0, 2048, CLMEM(debugBuffer), 0, NULL, NULL);
		}

		// Wait for async commands to finish
		CL(Finish, cmd);
		if (debugBuffer) if (handleDebugOutput()) exit(1);

		restartFinishedPhotons(photonCount, totalThreadCount, &finishedPhotonCount, R_specular);

		uint32_t totalFinishedPhotonCount = finishedPhotonCount;
		MPI(Reduce, &finishedPhotonCount, &totalFinishedPhotonCount, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);

		// Print status
		if (rank == 0) out << '\r' << "Photons terminated: " << totalFinishedPhotonCount << "/" << targetPhotonCount << Log::flush;
	}
	if (rank == 0) out << '\n';
}


// Declaration of kernel function to be able to use it from CPU code
void mcml(
	float nAbove,
	float nBelow,
	Layer* layers,
	int layerCount,
	Boundary* boundaries,
	int size_r,
	int size_a,
	int size_z,
	float delta_r, 
	float delta_z,
	volatile Weight* R_ra,
	volatile Weight* T_ra,
	volatile Weight* A_rz,
	PhotonTracker* photonStates);


static void simulateOnCPU(SimulationStruct sim, Layer* layers, Boundary* boundaries, Weight* R, Weight* A, Weight* T,
uint32_t photonCount, uint32_t targetPhotonCount, int rank, float R_specular) {
	float nAbove = sim.layers[0].n; // first and last layer have only n initialized...
	float nBelow = sim.layers[sim.n_layers + 1].n; // ...and represent outer media
	float radialBinCentimeters = sim.det.dr;
	float depthBinCentimeters = sim.det.dz;
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	int depthBinCount = sim.det.nz;
	uint32_t finishedPhotonCount = 0;
	if (rank == 0) out << '\n';
	while (finishedPhotonCount < photonCount) {
		mcml(
			nAbove,
			nBelow,
			layers,
			sim.n_layers,
			boundaries,
			radialBinCount,
			angularBinCount,
			depthBinCount,
			radialBinCentimeters,
			depthBinCentimeters,
			R, T, A,
			(PhotonTracker*)CLMEM(stateBuffer));
		restartFinishedPhotons(photonCount, 1, &finishedPhotonCount, R_specular);
		uint32_t totalFinishedPhotonCount = finishedPhotonCount;
		MPI(Reduce, &finishedPhotonCount, &totalFinishedPhotonCount, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
		// Print status
		if (rank == 0) out << '\r' << "Photons terminated: " << totalFinishedPhotonCount << "/" << targetPhotonCount << Log::flush;
	}
	if (rank == 0) out << '\n';
}


static void downloadOutputArrays(SimulationStruct sim, cl_command_queue cmd, cl_mem R, cl_mem A, cl_mem T) {
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	int depthBinCount = sim.det.nz;
	size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(Weight);
	size_t transmissionBufferSize = reflectanceBufferSize;
	size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(Weight);
	CL(EnqueueReadBuffer, cmd, R, CL_FALSE, 0, reflectanceBufferSize, CLMEM(R), 0, NULL, NULL);
	CL(EnqueueReadBuffer, cmd, T, CL_FALSE, 0, transmissionBufferSize, CLMEM(T), 0, NULL, NULL);
	CL(EnqueueReadBuffer, cmd, A, CL_FALSE, 0, absorptionBufferSize, CLMEM(A), 0, NULL, NULL);
	CL(Finish, cmd);
}


static void reduceOutputArrays(SimulationStruct sim, cl_mem R, cl_mem A, cl_mem T, Weight* outR, Weight* outA, Weight* outT) {
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	int depthBinCount = sim.det.nz;
	// Sum RAT buffers from all processes
	MPI(Reduce, CLMEM(R), outR, radialBinCount*angularBinCount, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI(Reduce, CLMEM(T), outT, radialBinCount*angularBinCount, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI(Reduce, CLMEM(A), outA, radialBinCount*depthBinCount, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
}


static void writeMCOFile(SimulationStruct sim, Weight* R_ra, Weight* A_rz, Weight* T_ra, cl_event kernelEvent) {
	uint64_t timeStart = 0, timeEnd = 0;
	CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(uint64_t), &timeStart, NULL);
	CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(uint64_t), &timeEnd, NULL);
	out << "Last Kerneltime=" << (timeEnd - timeStart) << "ns=" << (timeEnd - timeStart) / 1000000.0f << "ms\n";

	Write_Simulation_Results(A_rz, T_ra, R_ra, &sim, timeEnd - timeStart);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// You saw above the private implementation. Now follow the public interface routines.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* mcmlOptions, int rank) {
	// Read input file with process 0
	bool ignoreA = strstr(kernelOptions, "-D IGNORE_A") != NULL;
	if (rank == 0) readMCIFile(mcmlOptions, ignoreA, true, &simCount);

	#ifndef NO_GPU
	broadcastInputData(rank);
	#endif

	// Allocate per-simulation arrays
	layersPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	boundariesPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	reflectancePerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	transmissionPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	absorptionPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));

	// Photon states buffer
	stateBuffer = CLMALLOC_OUTPUT(totalThreadCount, PhotonTracker);

	// Allocate memory for each simulation
	for (int simIndex = 0; simIndex < simCount; simIndex++) {
		// ensure no bins can overflow
		assert(simulations[simIndex].number_of_photons <= 0xFFFFFFFFu);
		// RAT buffer mem
		allocOutputArrays(simulations[simIndex], simIndex);
		// layer and boundary mem
		setupInputArrays(simulations[simIndex], simIndex);
	}

	// Debug buffer
	int debugMode = strstr(kernelOptions, "-D DEBUG") != NULL ? 1 : 0;
	if (debugMode) debugBuffer = CLMALLOC_OUTPUT(2048, char);
}


void runCLKernel(cl_command_queue cmdQueue, cl_kernel kernel, size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		uploadInputArrays(cmdQueue, simulations[simIndex], layersPerSimulation[simIndex], boundariesPerSimulation[simIndex]);

		initAndUploadOutputArrays(cmdQueue, simulations[simIndex], reflectancePerSimulation[simIndex], absorptionPerSimulation[simIndex], transmissionPerSimulation[simIndex]);

		setKernelArguments(kernel, simulations[simIndex], layersPerSimulation[simIndex], boundariesPerSimulation[simIndex],
			reflectancePerSimulation[simIndex], absorptionPerSimulation[simIndex], transmissionPerSimulation[simIndex]);

		// Compute percentage of light that photons lose before entering the first layer
		float R_specular = computeFresnelReflectance(simulations[simIndex].layers[0].n, simulations[simIndex].layers[1].n);

		initPhotonStates(totalThreadCount, R_specular);

		// Distribute photons on processes
		uint32_t targetPhotonCount = simulations[simIndex].number_of_photons;
		uint32_t processPhotonCount = targetPhotonCount / processCount;
		if (rank == 0) { // rank 0 takes the remainder on top
			processPhotonCount += (targetPhotonCount % processCount);
		}
		cl_event kernelEvent;
		#ifndef NO_GPU
		simulate(cmdQueue, kernel, totalThreadCount, simdThreadCount, processPhotonCount, targetPhotonCount, rank, R_specular, kernelEvent);
		#else
		simulateOnCPU(simulations[simIndex],
			(Layer*)CLMEM(layersPerSimulation[simIndex]),
			(Boundary*)CLMEM(boundariesPerSimulation[simIndex]),
			(Weight*)CLMEM(reflectancePerSimulation[simIndex]),
			(Weight*)CLMEM(absorptionPerSimulation[simIndex]),
			(Weight*)CLMEM(transmissionPerSimulation[simIndex]),
			processPhotonCount, targetPhotonCount, rank, R_specular);
		#endif

		downloadOutputArrays(simulations[simIndex], cmdQueue, reflectancePerSimulation[simIndex], absorptionPerSimulation[simIndex], transmissionPerSimulation[simIndex]);

		#ifndef NO_GPU
		int radialBinCount = simulations[simIndex].det.nr;
		int angularBinCount = simulations[simIndex].det.na;
		int depthBinCount = simulations[simIndex].det.nz;
		size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(Weight);
		size_t transmissionBufferSize = reflectanceBufferSize;
		size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(Weight);
		Weight* totalReflectance = (Weight*)malloc(reflectanceBufferSize);
		Weight* totalTransmittance = (Weight*)malloc(transmissionBufferSize);
		Weight* totalAbsorption = (Weight*)malloc(absorptionBufferSize);
		reduceOutputArrays(simulations[simIndex], reflectancePerSimulation[simIndex], absorptionPerSimulation[simIndex], transmissionPerSimulation[simIndex],
			totalReflectance, totalTransmittance, totalAbsorption);
		#else
		Weight* totalReflectance = (Weight*)CLMEM(reflectancePerSimulation[simIndex]);
		Weight* totalTransmittance = (Weight*)CLMEM(absorptionPerSimulation[simIndex]);
		Weight* totalAbsorption = (Weight*)CLMEM(transmissionPerSimulation[simIndex]);
		#endif

		if (rank == 0) writeMCOFile(simulations[simIndex], totalReflectance, totalAbsorption, totalTransmittance, kernelEvent);

		#ifndef NO_GPU
		free(totalReflectance);
		free(totalTransmittance);
		free(totalAbsorption);
		#endif

		CL(ReleaseEvent, kernelEvent);
	}
	freeResources();
}


const char* getCLKernelName() {
	return "mcml";
}