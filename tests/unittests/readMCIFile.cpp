#include <iostream>
#include <assert.h>

#include "CUDAMCMLio.h" // SimulationStruct
#include "CUDAMCMLio.c" // read_simulation_data, Write_Simulation_Results

static int simCount;
static SimulationStruct* simulations;


/**
* Read data from mci file into array of simulation structs.
*/
static void readMCIFile(char* name, bool ignoreA, int* outSimCount) {

	*outSimCount = read_simulation_data(name, &simulations, ignoreA ? 1 : 0);

	assert(*outSimCount > 0);
}


static void dumpStructs() {
	SimulationStruct sim = simulations[0];

	for (int i = 0; i < sim.n_layers+1; i++) {
		std::cout << "BoundaryStruct:\n"
			<< "isHeightfield = " << sim.boundaries[i].isHeightfield << '\n';
		if (sim.boundaries[i].isHeightfield) {
			std::cout << "n=" << sim.boundaries[i].n << "\n";
			std::cout << "heights = ";
			for (int j = 0; j < sim.boundaries[i].n; j++) {
				std::cout << sim.boundaries[i].heights[j] << " ";
			}
			std::cout << "\nspacings = ";
			for (int j = 0; j < sim.boundaries[i].n; j++) {
				std::cout << sim.boundaries[i].spacings[j] << " ";
			}
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	for (int i = 0; i < sim.n_layers+2; i++) {
		std::cout << "LayerStruct:\n";
		std::cout << "n = " << sim.layers[i].n << std::endl;
		if (i == 0 || i == sim.n_layers+1) {
			continue;
		}
		std::cout << "z_min = " << sim.layers[i].z_min << std::endl;
		std::cout << "z_max = " << sim.layers[i].z_max << std::endl;
		std::cout << "mutr = " << sim.layers[i].mutr << std::endl;
		std::cout << "mua = " << sim.layers[i].mua << std::endl;
		std::cout << "g = " << sim.layers[i].g << std::endl;
	}
}



int main(int nargs, char** args) {

	char* file = args[1];

	std::cout << "readMCIFile( " << file << ", false, ...)\n\n";

	readMCIFile(file, false, &simCount);
	dumpStructs();

	return 0;
}