/*********************************************************************************
*
* This code sets up MCML and distributes work on processes and threads.
*
* It triggers a compute kernel until all photons are processed.
*
* It takes a single .mci-file with multiple simulations.
* It produces .mco-files for each simulation.
*
*********************************************************************************/

#include <string.h> // strstr, memcpy
#include <stdint.h> // uint32_t, uint64_t
#include <stdlib.h> // malloc, free
#include <assert.h> // assert

#include <time.h> // clock, CLOCKS_PER_SEC
#include <chrono> // std::chrono::high_resolution_clock::time_point, std::chrono::duration

// Memory management
#include "clmem.h" // CLMALLOC_INPUT, CLMALLOC_OUTPUT, CLMEM, CLMEM_ACCESS_ARRAY, CLMEM_ACCESS_AOS

// Core structs
#include "PhotonTracker.h" // struct PhotonTracker
#include "Layer.h" // struct Layer
#ifdef NO_GPU
typedef struct { float x,y,z; } float3;
#endif
#include "Boundary.h" // struct Boundary

// RNG functions
#include "randomlib.h" // wang_hash, rand_xorshift

// Logging class
#include "Log.h" // out

// Error catching macros
#define DEBUG // debug mode, i.e. enabled checks, always better for now
#include "clcheck.h" // CL macro
#include "mpicheck.h" // MPI macro

// MCML file io
//TODO Link instead of include
#include "CUDAMCMLio.h" // SimulationStruct
#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

// Photon weight type
typedef uint64_t Weight;
#define MPI_WEIGHT_T MPI_UINT64_T


/**
* MPI message type
*/
enum { GO, CONTINUE, FINISH };


/**
* Parallelization strategies
* UPFRONT_WORK_SPLIT is coarse-grained
* COORDINATOR_WORKER is fine-grained
*/
enum ParallelizationStrategy { UPFRONT_WORK_SPLIT, COORDINATOR_WORKER };

// Holds current parallelization strategy
const ParallelizationStrategy CURRENT_PAR = ParallelizationStrategy::COORDINATOR_WORKER;


// Holds a description of simulations, equivalent to mci file
static int simCount = 0;
static SimulationStruct* simulations = 0;

// Per-simulation buffers
static cl_mem* layersPerSimulation = 0;
static cl_mem* boundariesPerSimulation = 0;
static cl_mem* heightsPerSimulation = 0;
static cl_mem* spacingsPerSimulation = 0;
static cl_mem* reflectancePerSimulation = 0;
static cl_mem* transmissionPerSimulation = 0;
static cl_mem* absorptionPerSimulation = 0;

// Holds photon state for every GPU thread
static cl_mem stateBuffer = 0;

// Ignore explicit heightfield boundaries, even if present in input file?
bool ignoreHeightfields = false; // set via kernel define

// Holds photon ages in GPU runs without restart (for debug)
static uint32_t* ages;

// Holds random seed for each photon in its current GPU run (for debug)
static uint32_t* seeds;

// Debug buffer can be used to output error messages from GPU
// (only allocated if kernel is compiled with -D DEBUG flag)
static cl_mem debugBuffer = 0;


/**
* Get data for photon in its start configuration
*/
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


/**
* Check if explicitly defined heightfield boundaries overlap.
* If so, the boundaries are invalid and the program is terminated.
*/
static void checkBoundaries(BoundaryStruct* boundaries, int n) {
	//TODO
}


/**
* Read data from mci file into array of simulation structs.
*/
static void readMCIFile(char* name, bool ignoreA, int* outSimCount) {

	*outSimCount = read_simulation_data(name, &simulations, ignoreA ? 1 : 0);

	assert(*outSimCount > 0);
}


//TODO write automatic parameter-combination- and test-script
//     (also need a way to visualize curves to detect noise, which means that photon count must be increased)


/**
* Create simulation struct description to enable transfer via MPI
*/
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


/**
* Create layer struct description to enable transfer via MPI
*/
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


/**
* Create boundary struct description to enable transfer via MPI
*/
static MPI_Datatype createMPIBoundaryStruct() {
	// Describing the following data layout:
	// uint32_t isHeightfield;
	// uint32_t n;
	// float* heights;
	// float* spacings;
	int blockLengths[2] = {2, 2*sizeof(void*)};
	int offsets[2] = {0, 2*sizeof(uint32_t)};
	MPI_Datatype types[2] = {MPI_UINT32_T, MPI_CHAR};
	MPI_Datatype boundaryStruct;
	MPI(Type_create_struct, 2, blockLengths, offsets, types, &boundaryStruct);
	MPI(Type_commit, &boundaryStruct);
	return boundaryStruct;
}


/**
* Transfer all the input structs via MPI to all cluster nodes.
* This includes simulation data with pointers to layer and boundary data.
*/
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

	// Broadcast layer struct
	MPI_Datatype mpiLayerStruct = createMPILayerStruct();
	int mpiLayerSize = 0;
	MPI(Type_size, mpiLayerStruct, &mpiLayerSize);
	assert(mpiLayerSize == sizeof(LayerStruct));
	for (int i = 0; i < simCount; i++) {
		if (rank != 0)
			simulations[i].layers = (LayerStruct*)malloc((simulations[i].n_layers + 2) * sizeof(LayerStruct));
		MPI(Bcast, simulations[i].layers, simulations[i].n_layers + 2, mpiLayerStruct, 0, MPI_COMM_WORLD);
	}

	// Broadcast boundary struct
	MPI_Datatype mpiBoundaryStruct = createMPIBoundaryStruct();
	int mpiBoundarySize = 0;
	MPI(Type_size, mpiBoundaryStruct, &mpiBoundarySize);
	assert(mpiBoundarySize == sizeof(BoundaryStruct));
	for (int i = 0; i < simCount; i++) {
		if (rank != 0)
			simulations[i].boundaries = (BoundaryStruct*)malloc((simulations[i].n_layers + 1) * sizeof(BoundaryStruct));
		MPI(Bcast, simulations[i].boundaries, simulations[i].n_layers + 1, mpiBoundaryStruct, 0, MPI_COMM_WORLD);
	}

	// Broadcast heightfield data for boundaries that are heightfields
	for (int i = 0; i < simCount; i++) {
		uint32_t n_boundaries = simulations[i].n_layers + 1;
		for (int j = 0; j < n_boundaries; j++) {
			if (simulations[i].boundaries[j].isHeightfield) {
				uint32_t n_samples = simulations[i].boundaries[j].n;
				if (rank != 0) {
					simulations[i].boundaries[j].heights = (float*)malloc(n_samples * sizeof(float));
					simulations[i].boundaries[j].spacings = (float*)malloc(n_samples * sizeof(float));
				}
				MPI(Bcast, simulations[i].boundaries[j].heights, n_samples, MPI_FLOAT, 0, MPI_COMM_WORLD);
				MPI(Bcast, simulations[i].boundaries[j].spacings, n_samples, MPI_FLOAT, 0, MPI_COMM_WORLD);
			}
		}
	}
}



/**
* Restructure layer and boundary data to be ready for GPU consumption.
* This takes data from given SimulationStruct to fill the separate buffers.
* Buffer handles are held by this module in {layers|boundaries}PerSimulation at the given index.
*/
static void setupInputArrays(SimulationStruct sim, int simIndex) {

	// Alloc space for layers and boundaries
	int layerCount = sim.n_layers;
	layersPerSimulation[simIndex] = CLMALLOC_INPUT(layerCount, Layer);
	boundariesPerSimulation[simIndex] = CLMALLOC_INPUT(layerCount+1, Boundary);
	
	unsigned int totalSampleCount = 0;
	for (int i = 0; i < layerCount+1; i++)
		if (sim.boundaries[i].isHeightfield && !ignoreHeightfields)
			totalSampleCount += sim.boundaries[i].n;

	heightsPerSimulation[simIndex] = CLMALLOC_INPUT(totalSampleCount, float);
	spacingsPerSimulation[simIndex] = CLMALLOC_INPUT(totalSampleCount, float);

	unsigned int sampleCount = 0;

	// Transform data into format that we want to use on GPU
	for (int j = 1; j <= layerCount; j++) {

		// Boundary
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, isHeightfield) = 
			sim.boundaries[j-1].isHeightfield && !ignoreHeightfields;
		if (CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, isHeightfield)) {
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.center.x) = 0.0f;
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.center.y) = 0.0f;
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.center.z) = sim.layers[j].z_min;

			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.i_heights) = sampleCount;
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.n_heights) = sim.boundaries[j-1].n;
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.i_spacings) = sampleCount;
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, heightfield.n_spacings) = sim.boundaries[j-1].n;
			memcpy(
				&CLMEM_ACCESS_ARRAY(CLMEM(heightsPerSimulation[simIndex]), float, sampleCount),
				sim.boundaries[j-1].heights,
				sim.boundaries[j-1].n * sizeof(float));
			memcpy(
				&CLMEM_ACCESS_ARRAY(CLMEM(spacingsPerSimulation[simIndex]), float, sampleCount),
				sim.boundaries[j-1].spacings,
				sim.boundaries[j-1].n * sizeof(float));

			sampleCount += sim.boundaries[j-1].n;
		} else {
			CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, j-1, z) = sim.layers[j].z_min;
		}

		// Layer
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, absorbCoeff) = sim.layers[j].mua;
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, scatterCoeff) = 
			1.0f / sim.layers[j].mutr - sim.layers[j].mua;
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, g) = sim.layers[j].g;
		CLMEM_ACCESS_AOS(CLMEM(layersPerSimulation[simIndex]), Layer, j-1, n) = sim.layers[j].n;
	}

	// One more boundary
	CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, isHeightfield) = 
		sim.boundaries[layerCount].isHeightfield && !ignoreHeightfields;
	if (CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, isHeightfield)) {
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.center.x) = 0.0f;
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.center.y) = 0.0f;
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.center.z) = sim.layers[layerCount].z_max;

		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.i_heights) = sampleCount;
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.n_heights) = sim.boundaries[layerCount].n;
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.i_spacings) = sampleCount;
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, heightfield.n_spacings) = sim.boundaries[layerCount].n;

		memcpy(
			&CLMEM_ACCESS_ARRAY(CLMEM(heightsPerSimulation[simIndex]), float, sampleCount),
			sim.boundaries[layerCount].heights,
			sim.boundaries[layerCount].n * sizeof(float));
		memcpy(
			&CLMEM_ACCESS_ARRAY(CLMEM(spacingsPerSimulation[simIndex]), float, sampleCount),
			sim.boundaries[layerCount].spacings,
			sim.boundaries[layerCount].n * sizeof(float));
		// Note: To avoid bugs in this kind of copy operation, there should be size-checked functions in clmem.h
	} else {
		CLMEM_ACCESS_AOS(CLMEM(boundariesPerSimulation[simIndex]), Boundary, layerCount, z) = sim.layers[layerCount].z_max;
	}
}


/**
* Alloc mem for ouput data.
* Array dimensions are taken from given SimulationStruct.
* Mem handles are stored in {reflectance|absorption|transmittance}PerSimulation at given index.
*/
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


/**
* Transfer layer and boundary data to GPU
*/
static void uploadInputArrays(cl_command_queue cmd, SimulationStruct sim, 
	cl_mem layers, cl_mem boundaries, cl_mem heights, cl_mem spacings) {

	// Do a blocking write to be safe (but maybe slow)
	//TODO can we use async buffer transfers?

	CL(EnqueueWriteBuffer, cmd, layers, CL_TRUE, 0, sim.n_layers*sizeof(Layer), CLMEM(layers), 0, NULL, NULL);
	CL(EnqueueWriteBuffer, cmd, boundaries, CL_TRUE, 0, (sim.n_layers+1)*sizeof(Boundary), CLMEM(boundaries), 0, NULL, NULL);

	unsigned int totalSampleCount = 0;
	for (int i = 0; i < sim.n_layers+1; i++)
		if (sim.boundaries[i].isHeightfield && !ignoreHeightfields)
			totalSampleCount += sim.boundaries[i].n;

	CL(EnqueueWriteBuffer, cmd, heights, CL_TRUE, 0, totalSampleCount*sizeof(float), CLMEM(heights), 0, NULL, NULL);
	CL(EnqueueWriteBuffer, cmd, spacings, CL_TRUE, 0, totalSampleCount*sizeof(float), CLMEM(spacings), 0, NULL, NULL);
}


/**
* Transfer zero-initialized output arrays to GPU
*/
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


/**
* Get fresnel reflectance at interface of two media with refractive indices n0 and n1
*/
static float computeFresnelReflectance(float n0, float n1) {
	// Reflectance (specular):
	// percentage of light leaving at surface without any interaction
	// using MCML's Fesnel approximation similar to Schlick's approximation (no incident angle, no polarization)
	// Q: why are the ^2 different than in Schlicks approximation?
	float nDiff = n0 - n1;
	float nSum = n0 + n1;
	return (nDiff * nDiff) / (nSum * nSum);
}


/**
* Init tracking data for totalThreadCount photons.
* This includes weight loss through fresnel reflection at first layer and
* seeding the RNG.
*/
static void initPhotonStates(size_t totalThreadCount, float R_specular) {
	for (unsigned int i = 0; i < totalThreadCount; i++) {
		PhotonTracker newState = createNewPhotonTracker();
		newState.weight -= R_specular;
		newState.rngState = wang_hash(i); // Seed RNG

		seeds[i] = newState.rngState;

		CLMEM_ACCESS_ARRAY(CLMEM(stateBuffer), PhotonTracker, i) = newState;
	}
}


/**
* Restart finished photons while ensuring not to start more than the photon count assigned to this MPI process.
* This function interates over all totalThreadCount photons (sequentially) to check for zero weight.
*/
static void restartFinishedPhotons(uint32_t processPhotonCount, size_t totalThreadCount, uint32_t* outFinishCount, float R_specular) {
	for (unsigned int i = 0; i < totalThreadCount; i++) {

		if (!CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead)) {

			ages[i]++; // this photon got one gpu round older

			if (CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) == 0) { // has it died?

				(*outFinishCount)++;

				// Launch new photon only if next round cannot overachieve
				if ((*outFinishCount) + totalThreadCount <= processPhotonCount) {

					ages[i] = 0;

					PhotonTracker newState = createNewPhotonTracker();
					newState.weight -= R_specular;
					newState.rngState = wang_hash((*outFinishCount) + totalThreadCount);

					seeds[i] = newState.rngState;

					CLMEM_ACCESS_ARRAY(CLMEM(stateBuffer), PhotonTracker, i) = newState;

				} else {
					CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead) = 1;
				}
			}
		}
	}
}



// ==========================================================================================
//TODO fuse the following functions into one "photon state analyzing and restart"-step

static uint32_t countFinishedPhotons(cl_command_queue cmd, size_t totalThreadCount) {

	// Download photon states
	CL(EnqueueReadBuffer, cmd, stateBuffer, CL_FALSE, 0,
		totalThreadCount * sizeof(PhotonTracker), CLMEM(stateBuffer), 0, NULL, NULL);

	uint32_t count = 0;

	//TODO omp parallel for
	for (unsigned int i = 0; i < totalThreadCount; i++) {
		if (!CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead)) {
			if (CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) == 0) {
				count++;
			}
		}
	}

	return count;
}

static void respawnFinishedPhotons(cl_command_queue cmd, size_t totalThreadCount, float R_specular) {

	//TODO omp parallel for
	for (unsigned int i = 0; i < totalThreadCount; i++) {
		if (!CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead)) {
			if (CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) == 0) {

					PhotonTracker newState = createNewPhotonTracker();
					newState.weight -= R_specular;

					// Reuse last rng state
					newState.rngState = CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, rngState);

					CLMEM_ACCESS_ARRAY(CLMEM(stateBuffer), PhotonTracker, i) = newState;
			}
		}
	}

	// Upload photon states
	//TODO since buffer updates are sparse, map could be faster than write in whole
	CL(EnqueueWriteBuffer, cmd, stateBuffer, CL_FALSE, 0,
		totalThreadCount * sizeof(PhotonTracker), CLMEM(stateBuffer), 0, NULL, NULL);
}

static void spawnExactPhotonCount(cl_command_queue cmd, size_t totalThreadCount, float R_specular, uint32_t count) {

	//TODO omp parallel for
	for (unsigned int i = 0; i < totalThreadCount; i++) {

		if (i >= count) {
			CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) = 0;
			CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead) = 1;

		} else {

			PhotonTracker newState = createNewPhotonTracker();
			newState.weight -= R_specular;

			// Reuse last rng state
			newState.rngState = CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, rngState);

			CLMEM_ACCESS_ARRAY(CLMEM(stateBuffer), PhotonTracker, i) = newState;
		}
	}

	// Upload photon states
	//TODO since buffer updates are sparse, map could be faster than write in whole
	CL(EnqueueWriteBuffer, cmd, stateBuffer, CL_FALSE, 0,
		totalThreadCount * sizeof(PhotonTracker), CLMEM(stateBuffer), 0, NULL, NULL);
}

static void killFinishedPhotons(cl_command_queue cmd, size_t totalThreadCount) {

	//TODO omp parallel for
	for (unsigned int i = 0; i < totalThreadCount; i++) {
		if (CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) == 0) {
			CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, isDead) = 1;
		}
	}

	// Upload photon states
	//TODO since buffer updates are sparse, map could be faster than write in whole
	CL(EnqueueWriteBuffer, cmd, stateBuffer, CL_FALSE, 0,
		totalThreadCount * sizeof(PhotonTracker), CLMEM(stateBuffer), 0, NULL, NULL);
}

static uint32_t countActivePhotons(cl_command_queue cmd, size_t totalThreadCount) {

	// Download photon states
	CL(EnqueueReadBuffer, cmd, stateBuffer, CL_FALSE, 0,
		totalThreadCount * sizeof(PhotonTracker), CLMEM(stateBuffer), 0, NULL, NULL);

	uint32_t count = 0;

	//TODO omp parallel for
	for (unsigned int i = 0; i < totalThreadCount; i++) {
		if (CLMEM_ACCESS_AOS(CLMEM(stateBuffer), PhotonTracker, i, weight) == 0) {
			count++;
		}
	}

	return totalThreadCount - count;
}


// ==========================================================================================


/**
* Set arguments for CL kernel function
*/
static void setKernelArguments(cl_kernel kernel, SimulationStruct sim, 
	cl_mem layers, cl_mem boundaries, cl_mem heights, cl_mem spacings, cl_mem R, cl_mem A, cl_mem T) {

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
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &heights); // heightfield samples for all boundaries
	CL(SetKernelArg, kernel, argCount++, sizeof(cl_mem), &spacings); // spacings between samples
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




// ==========================================================================================

// Fine grained parallelization approach (MPI coordinator/worker pattern)

/**
*
*/
static void runMPICoordinator(SimulationStruct* sim, cl_command_queue cmdQueue, cl_kernel kernel, 
	size_t totalThreadCount, size_t simdThreadCount, float R_specular, int processCount) {

	// Number of unstarted photons
	uint32_t L = sim->number_of_photons;

	// Number of workers finished
	int finishCount = 0;

	// Tell workers to start
	// (This is a necessary synchronization step, because 
	// it turned out that when workers directly start by asking for work
	// the corresponding send call hangs indefinitely.
	// Usually this happens when the other side is not in receiving state.)
	for (int i = 1; i < processCount; i++)
		MPI(Send, 0, 0, MPI_INT, i, GO, MPI_COMM_WORLD);
	
	while (finishCount < (processCount-1)) { // not all workers finished?

		MPI_Status status;

		// Requested photons
		uint32_t K;

		// Wait for request
		MPI(Recv, &K, 1, MPI_UINT32_T, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		// Respond
		if (L >= K) {
			MPI(Send, 0, 0, MPI_INT, status.MPI_SOURCE, CONTINUE, MPI_COMM_WORLD);
			L -= K;
		} else {
			MPI(Send, 0, 0, MPI_INT, status.MPI_SOURCE, FINISH, MPI_COMM_WORLD);
			finishCount++;
		}
	}

	if (L > 0) {

		spawnExactPhotonCount(cmdQueue, totalThreadCount, R_specular, L);

		while (L > 0) {
			killFinishedPhotons(cmdQueue, totalThreadCount);
			CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, 0);
			CL(Finish, cmdQueue);
			L -= countFinishedPhotons(cmdQueue, totalThreadCount);
		}
	}
}

/**
*
*/
static void runMPIWorker(SimulationStruct* sim, cl_command_queue cmdQueue, cl_kernel kernel, 
	size_t totalThreadCount, size_t simdThreadCount, float R_specular, int rank) {

	// Communication threshold
	// 0 means that if only one photon finished a request is made to restart it
	// higher number means buying lower communication overhead with higher number of inactive threads
	const uint32_t Z = 0;

	// Current finished photons, i.e. inactive threads
	uint32_t K = (uint32_t)totalThreadCount;

	MPI_Status status;

	// Wait for GO message
	MPI(Recv, 0, 0, MPI_INT, 0, GO, MPI_COMM_WORLD, &status);

	while (true) {

		if (K > Z) {

			// Request to respawn K photons
			MPI(Send, &K, 1, MPI_UINT32_T, 0, 0, MPI_COMM_WORLD);

			// Wait for response
			MPI(Recv, 0, 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if (status.MPI_TAG == CONTINUE) {

				// Respawn K photons
				respawnFinishedPhotons(cmdQueue, totalThreadCount, R_specular);
				K = 0;

			} else if (status.MPI_TAG == FINISH) {

				// Finish remaining photons without further requests
				// (more and more threads will become inactive to meet exact photon count,
				// which can be optimized by instant restart on GPU with atomic finish counter)

				bool isFinished = false;

				while (!isFinished) {
					killFinishedPhotons(cmdQueue, totalThreadCount);
					CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, 0);
					CL(Finish, cmdQueue);
					uint32_t tmp = countActivePhotons(cmdQueue, totalThreadCount);
					isFinished = (tmp == 0);
				}

				break;

			}
		}

		// Continue processing
		CL(EnqueueNDRangeKernel, cmdQueue, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, 0);
		CL(Finish, cmdQueue);
		K += countFinishedPhotons(cmdQueue, totalThreadCount);
	}
}




// ==========================================================================================

// Coarse grained parallelization approach: work was split up-front and 
// each process runs simulate() independently


/**
* Do simulation on GPU for all photons assigned to this MPI process
*/
static void simulate(cl_command_queue cmd, cl_kernel kernel, size_t totalThreadCount, size_t simdThreadCount, 
uint32_t photonCount, uint32_t targetPhotonCount, int rank, float R_specular) {
	assert(photonCount >= totalThreadCount); // this is a non-sensical low photon count

	size_t photonStateBufferSize = totalThreadCount * sizeof(PhotonTracker);

	uint32_t finishedPhotonCount = 0;
	uint32_t gpuRoundCounter = 0;

	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	cl_event kernelEvent;
	uint64_t totalKernelTime = 0; // nanoseconds
	cl_event transferEvent;
	uint64_t totalTransferTime = 0; // nanoseconds

	if (rank == 0) out << '\n';

	while (finishedPhotonCount < photonCount) {

		// Upload photon states
		//TODO since buffer updates are sparse, map could be faster than write in whole
		CL(EnqueueWriteBuffer, cmd, stateBuffer, CL_FALSE, 0,
			photonStateBufferSize, CLMEM(stateBuffer), 0, NULL, NULL);

		// Run a batch of photons
		CL(EnqueueNDRangeKernel, cmd, kernel, 1, NULL, &totalThreadCount, &simdThreadCount, 0, NULL, &kernelEvent);
		gpuRoundCounter++;

		// Download photon states
		CL(EnqueueReadBuffer, cmd, stateBuffer, CL_FALSE, 0, photonStateBufferSize, CLMEM(stateBuffer), 0, NULL, &transferEvent);
		if (debugBuffer) { // download debug buffer if in debug mode
			CL(EnqueueReadBuffer, cmd, debugBuffer, CL_FALSE, 0, 2048, CLMEM(debugBuffer), 0, NULL, NULL);
		}

		// Wait for async commands to finish
		CL(Finish, cmd);

		uint64_t timeStart, timeEnd;
		CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(uint64_t), &timeStart, NULL);
		CL(GetEventProfilingInfo, kernelEvent, CL_PROFILING_COMMAND_END, sizeof(uint64_t), &timeEnd, NULL);
		totalKernelTime += (timeEnd - timeStart);
		CL(GetEventProfilingInfo, transferEvent, CL_PROFILING_COMMAND_QUEUED, sizeof(uint64_t), &timeStart, NULL);
		CL(GetEventProfilingInfo, transferEvent, CL_PROFILING_COMMAND_END, sizeof(uint64_t), &timeEnd, NULL);
		totalTransferTime += (timeEnd - timeStart);

		restartFinishedPhotons(photonCount, totalThreadCount, &finishedPhotonCount, R_specular);

		// Synchronize to report total progress
		// (can introduce much waiting time, especially when done in each iteration)
		uint32_t totalFinishedPhotonCount = finishedPhotonCount;
		MPI(Reduce, &finishedPhotonCount, &totalFinishedPhotonCount, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0)
			out << '\r' << "Photons terminated: " << totalFinishedPhotonCount <<
				"/" << targetPhotonCount << 
				", Round " << gpuRoundCounter <<
				", Time elapsed = " << std::chrono::duration_cast<std::chrono::nanoseconds>(
					std::chrono::high_resolution_clock::now() - start).count() / 1000000.0 << "ms" <<
				", Kernel time = " << totalKernelTime/1000000.0 << "ms" <<
				", Transfer time = " << totalTransferTime/1000000.0 << "ms" << Log::flush;

	}
	if (rank == 0) out << '\n';

	//TODO compare perf against CUDAMCML, which does not wait for longest simulating thread,
	// but instead launches a new photon directly from a thread that would terminate,
	// causing additional tracking overhead.
}



// ==========================================================================================
// NO_GPU and single process: completely sequential and easier to debug
// Should also behave exactly like original mcml (same RNG sequence)


/**
* Declaration of kernel function to be able to use it from CPU code
*/
void mcml(
	float nAbove,
	float nBelow,
	Layer* layers,
	int layerCount,
	Boundary* boundaries,
	float* heights,
	float* spacings,
	int size_r,
	int size_a,
	int size_z,
	float delta_r, 
	float delta_z,
	volatile Weight* R_ra,
	volatile Weight* T_ra,
	volatile Weight* A_rz,
	PhotonTracker* photonStates);


/**
* Do simulation on CPU, for debug purposes
*/
static void simulateOnCPU(SimulationStruct sim,
	Layer* layers, Boundary* boundaries, float* heights, float* spacings, Weight* R, Weight* A, Weight* T,
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
			heights,
			spacings,
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





// ==========================================================================================




/**
* Transfer output data back from GPU (synchronous transfer)
*/
static void downloadOutputArrays(SimulationStruct sim, cl_command_queue cmd, cl_mem R, cl_mem A, cl_mem T) {
	int radialBinCount = sim.det.nr, angularBinCount = sim.det.na, depthBinCount = sim.det.nz;
	size_t reflectanceBufferSize = radialBinCount * angularBinCount * sizeof(Weight);
	size_t transmissionBufferSize = reflectanceBufferSize;
	size_t absorptionBufferSize = radialBinCount * depthBinCount * sizeof(Weight);

	CL(EnqueueReadBuffer, cmd, R, CL_FALSE, 0, reflectanceBufferSize, CLMEM(R), 0, NULL, NULL);
	CL(EnqueueReadBuffer, cmd, A, CL_FALSE, 0, absorptionBufferSize, CLMEM(A), 0, NULL, NULL);
	CL(EnqueueReadBuffer, cmd, T, CL_FALSE, 0, transmissionBufferSize, CLMEM(T), 0, NULL, NULL);
	if (debugBuffer) CL(EnqueueReadBuffer, cmd, debugBuffer, CL_FALSE, 0, 2048, CLMEM(debugBuffer), 0, NULL, NULL);
	CL(Finish, cmd);
}


/**
* Sum up the weights in the output data over all MPI processes
*/
static void reduceOutputArrays(SimulationStruct sim, cl_mem R, cl_mem A, cl_mem T, Weight* outR, Weight* outA, Weight* outT) {
	int radialBinCount = sim.det.nr;
	int angularBinCount = sim.det.na;
	int depthBinCount = sim.det.nz;
	// Sum RAT buffers from all processes
	MPI(Reduce, CLMEM(R), outR, radialBinCount*angularBinCount, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI(Reduce, CLMEM(T), outT, radialBinCount*angularBinCount, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI(Reduce, CLMEM(A), outA, radialBinCount*depthBinCount, MPI_WEIGHT_T, MPI_SUM, 0, MPI_COMM_WORLD);
}


/**
* Write mco file
*/
static void writeMCOFile(SimulationStruct sim, Weight* R_ra, Weight* A_rz, Weight* T_ra) {

	uint64_t timeStart = 0, timeEnd = 0; //TODO measure total elapsed time

	assert(Write_Simulation_Results(A_rz, T_ra, R_ra, &sim, timeEnd - timeStart));
	
	out << "Output file written: " << sim.outp_filename << "\n";
}


/**
* Release memory of everything
*/
static void freeResources() {
	for (int i = 0; i < simCount; i++) {
		free(simulations[i].layers);
		for (int j = 0; j < simulations[i].n_layers+1; j++) {
			if (simulations[i].boundaries[j].isHeightfield) {
				free(simulations[i].boundaries[j].heights);
				free(simulations[i].boundaries[j].spacings);
			}
		}
		free(simulations[i].boundaries);
	}
	free(absorptionPerSimulation);
	free(transmissionPerSimulation);
	free(reflectancePerSimulation);
	free(spacingsPerSimulation);
	free(heightsPerSimulation);
	free(boundariesPerSimulation);
	free(layersPerSimulation);
	free(simulations);
	free(ages);
	free(seeds);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// You saw above the private implementation. Now follows the public interface.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
* Setup all the stuff needed for CL+MPI implementation of MCML
*/
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* mcmlOptions, int rank) {

	bool ignoreA = strstr(kernelOptions, "-D IGNORE_A") != NULL;
	if (rank == 0) readMCIFile(mcmlOptions, ignoreA, &simCount);

	ignoreHeightfields = strstr(kernelOptions, "-D IGNORE_HEIGHTFIELDS") != NULL;
	if (!ignoreHeightfields)
		for (int i = 0; i < simCount; i++)
			checkBoundaries(simulations[i].boundaries, simulations[i].n_layers + 1);

	#ifndef NO_GPU
	{
		clock_t start = clock();
		broadcastInputData(rank);
		clock_t end = clock();
		out << "Time spent in MPI broadcast: " << (((double) (end - start)) / CLOCKS_PER_SEC)*1000 << "ms\n";
	}
	#endif

	// Allocate per-simulation arrays
	layersPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	boundariesPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	heightsPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	spacingsPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	reflectancePerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	transmissionPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	absorptionPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));

	// Photon states buffer
	stateBuffer = CLMALLOC_OUTPUT(totalThreadCount, PhotonTracker);

	ages = (uint32_t*)malloc(totalThreadCount * sizeof(uint32_t));
	seeds = (uint32_t*)malloc(totalThreadCount * sizeof(uint32_t));
	for (unsigned int i = 0; i < totalThreadCount; i++) {
		ages[i] = seeds[i] = 0;
	}

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
	if (debugMode) {
		out << "USING KERNEL IN DEBUG MODE\n";
		debugBuffer = CLMALLOC_OUTPUT(2048, char);
		for (unsigned int i = 0; i < 2048; i++)
			CLMEM_ACCESS_ARRAY(CLMEM(debugBuffer), char, i) = 0;
	}
}


/**
* Run all the simulations as specified in mci file
*/
void runCLKernel(cl_command_queue cmdQueue, cl_kernel kernel, size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank) {
	for (int simIndex = 0; simIndex < simCount; simIndex++) {

		uploadInputArrays(cmdQueue, simulations[simIndex], 
			layersPerSimulation[simIndex], 
			boundariesPerSimulation[simIndex],
			heightsPerSimulation[simIndex],
			spacingsPerSimulation[simIndex]);

		initAndUploadOutputArrays(cmdQueue, simulations[simIndex], 
			reflectancePerSimulation[simIndex], 
			absorptionPerSimulation[simIndex], 
			transmissionPerSimulation[simIndex]);

		setKernelArguments(kernel, simulations[simIndex], 
			layersPerSimulation[simIndex], 
			boundariesPerSimulation[simIndex],
			heightsPerSimulation[simIndex],
			spacingsPerSimulation[simIndex],
			reflectancePerSimulation[simIndex], 
			absorptionPerSimulation[simIndex], 
			transmissionPerSimulation[simIndex]);

		// Compute percentage of light that photons lose before entering the first layer
		float R_specular = computeFresnelReflectance(
			simulations[simIndex].layers[0].n, 
			simulations[simIndex].layers[1].n);

		// Initialize necessary start parameters
		initPhotonStates(totalThreadCount, R_specular);


		// Simulate

		#ifndef NO_GPU

			if (processCount > 1 && CURRENT_PAR == ParallelizationStrategy::COORDINATOR_WORKER) {

				if (rank == 0) {

					runMPICoordinator(&simulations[simIndex], cmdQueue, kernel, totalThreadCount, simdThreadCount, R_specular, processCount);

				} else {

					runMPIWorker(&simulations[simIndex], cmdQueue, kernel, totalThreadCount, simdThreadCount, R_specular, rank);
				}
			} else {

				uint32_t targetPhotonCount = simulations[simIndex].number_of_photons;
				uint32_t processPhotonCount = targetPhotonCount / processCount;
				if (rank == 0) { // rank 0 takes the remainder on top
					processPhotonCount += (targetPhotonCount % processCount);
				}
				simulate(cmdQueue, kernel, totalThreadCount, simdThreadCount, processPhotonCount, targetPhotonCount, rank, R_specular);

			}

		#else

			simulateOnCPU(simulations[simIndex],
				(Layer*)CLMEM(layersPerSimulation[simIndex]),
				(Boundary*)CLMEM(boundariesPerSimulation[simIndex]),
				(float*)CLMEM(heightsPerSimulation[simIndex]),
				(float*)CLMEM(spacingsPerSimulation[simIndex]),
				(Weight*)CLMEM(reflectancePerSimulation[simIndex]),
				(Weight*)CLMEM(absorptionPerSimulation[simIndex]),
				(Weight*)CLMEM(transmissionPerSimulation[simIndex]),
				simulations[simIndex].number_of_photons, simulations[simIndex].number_of_photons, rank, R_specular);
		
		#endif



		downloadOutputArrays(simulations[simIndex], cmdQueue, 
			reflectancePerSimulation[simIndex], absorptionPerSimulation[simIndex], transmissionPerSimulation[simIndex]);



		// Handle data that kernel wrote to debug buffer (kernel DEBUG define required)
		if (debugBuffer) {
			// Uncomment to print finished photon counter
			// uint32_t sum;
			// MPI(Reduce, CLMEM(debugBuffer), &sum, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);
			// if (rank == 0) out << "\nDEBUG PHOTON COUNT = " << sum << "\n\n";

			// Uncomment to print floats from debug buffer
			// for (int i = 0; i < 60; i++) {
			// 	out << ((float*)CLMEM(debugBuffer))[i] << ", ";
			// }
			// out << '\n';

			// Uncomment ot print parameters of first two boundaries
			// out << "isHeightfield = " << ((uint32_t*)CLMEM(debugBuffer))[0] << ", " <<((uint32_t*)CLMEM(debugBuffer))[1] << "\n";
			// out << "z = " << ((float*)CLMEM(debugBuffer))[2] << ", " << ((float*)CLMEM(debugBuffer))[3] << '\n';
			// out << "i_heights = " << ((uint32_t*)CLMEM(debugBuffer))[4] << ", " << ((uint32_t*)CLMEM(debugBuffer))[5] << '\n';
			// out << "n_heights = " << ((uint32_t*)CLMEM(debugBuffer))[6] << ", " << ((uint32_t*)CLMEM(debugBuffer))[7] << '\n';
			// out << "i_spacings = " << ((uint32_t*)CLMEM(debugBuffer))[8] << ", " << ((uint32_t*)CLMEM(debugBuffer))[9] << '\n';
			// out << "n_spacings = " << ((uint32_t*)CLMEM(debugBuffer))[10] << ", " << ((uint32_t*)CLMEM(debugBuffer))[11] << '\n';
			// out << "sizeof(Boundary) = " << ((uint32_t*)CLMEM(debugBuffer))[12] << '\n';
		}



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

		out << "Reduce... ";
		reduceOutputArrays(simulations[simIndex],
			reflectancePerSimulation[simIndex], absorptionPerSimulation[simIndex], transmissionPerSimulation[simIndex],
			totalReflectance, totalAbsorption, totalTransmittance);
		out << "done.\n";

		#else

		Weight* totalReflectance = (Weight*)CLMEM(reflectancePerSimulation[simIndex]);
		Weight* totalTransmittance = (Weight*)CLMEM(transmissionPerSimulation[simIndex]);
		Weight* totalAbsorption = (Weight*)CLMEM(absorptionPerSimulation[simIndex]);

		#endif



		if (rank == 0) {
			out << "Write to file...\n";
			writeMCOFile(simulations[simIndex], totalReflectance, totalAbsorption, totalTransmittance);
		}



		#ifndef NO_GPU
		free(totalReflectance);
		free(totalTransmittance);
		free(totalAbsorption);
		#endif

		// CL(ReleaseEvent, kernelEvent);
	}
	freeResources();
}


/**
* Get name of kernel function that does MCML photon transport
*/
const char* getCLKernelName() {
	return "mcml";
}