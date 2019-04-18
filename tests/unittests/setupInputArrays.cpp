#include <iostream>
#include <assert.h>

#include "clmem.h"

#include "CUDAMCMLio.h" // SimulationStruct
#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

#include "Layer.h" // struct Layer
typedef struct { float x,y,z; } float3;
#include "Boundary.h" // struct Boundary


// Data:

static int simCount;
static SimulationStruct* simulations;

static cl_mem* layersPerSimulation = 0;
static cl_mem* boundariesPerSimulation = 0;
static cl_mem* heightsPerSimulation = 0;
static cl_mem* spacingsPerSimulation = 0;



// Code:


/**
* Read data from mci file into array of simulation structs.
*/
static void readMCIFile(char* name, bool ignoreA, int* outSimCount) {

	*outSimCount = read_simulation_data(name, &simulations, ignoreA ? 1 : 0);

	assert(*outSimCount > 0);
}

/**
* Restructure layer and boundary data to be ready for GPU consumption.
* This takes data from given SimulationStruct to fill the separate buffers.
* Buffer handles are held by this module in {layers|boundaries}PerSimulation at the given index.
*/
static void setupInputArrays(SimulationStruct sim, int simIndex, bool ignoreHeightfields) {

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




static void dumpHeightfieldData(unsigned int boundaryCount, cl_mem boundaries, cl_mem heights, cl_mem spacings) {

	for (unsigned int i = 0; i < boundaryCount; i++) {
		
		RHeightfield field = CLMEM_ACCESS_AOS(CLMEM(boundaries), Boundary, i, heightfield);
		std::cout << std::endl << "Boundary " << i <<":" << std::endl;

		bool isHeightfield = (bool)CLMEM_ACCESS_AOS(CLMEM(boundaries), Boundary, i, isHeightfield);

		std::cout << "isHeightfield = " << isHeightfield << std::endl;

		if (!isHeightfield) {
			std::cout << "z = " << (float)CLMEM_ACCESS_AOS(CLMEM(boundaries), Boundary, i, z) << std::endl;
		}

		if (isHeightfield) {

			std::cout << "center = " << field.center.x << ", " << field.center.y << ", " << field.center.z << std::endl;
			std::cout << "i_heights = " << field.i_heights << std::endl;
			std::cout << "n_heights = " << field.n_heights << std::endl;
			std::cout << "i_spacings = " << field.i_spacings << std::endl;
			std::cout << "n_spacings = " << field.n_spacings << std::endl;

			std::cout << "heights = ";
			for (uint32_t i = field.i_heights; i < field.i_heights+field.n_heights; i++)
				std::cout << CLMEM_ACCESS_ARRAY(CLMEM(heights), float, i) << ", ";
			std::cout << std::endl;

			std::cout << "spacings = ";
			for (uint32_t i = field.i_spacings; i < field.i_spacings+field.n_spacings; i++)
				std::cout << CLMEM_ACCESS_ARRAY(CLMEM(spacings), float, i) << ", ";
			std::cout << std::endl;
		}
	}
}

int main(int nargs, char** args) {

	char* file = args[1];

	std::cout << "readMCIFile( " << file << ", false, ...)" << std::endl;

	readMCIFile(file, false, &simCount);

	// These arrays are to be filled
	layersPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	boundariesPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	heightsPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));
	spacingsPerSimulation = (cl_mem*)malloc(simCount * sizeof(cl_mem));

	for (int i = 0; i < simCount; i++) {

		std::cout << std::endl << "setupInputArrays(simulations["<<i<<"], "<<i<<", false)" << std::endl;

		setupInputArrays(simulations[i], i, false);

		dumpHeightfieldData(simulations[i].n_layers+1, boundariesPerSimulation[i], heightsPerSimulation[i], spacingsPerSimulation[i]);
	}

	return 0;
}