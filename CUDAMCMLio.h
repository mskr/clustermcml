/*	This file is part of CUDAMCML.

    CUDAMCML is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUDAMCML is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUDAMCML.  If not, see <http://www.gnu.org/licenses/>.*/

#ifndef CUDAMCML_H
#define CUDAMCML_H

#define PI 3.14159265359

#define STR_LEN 200 // should align to 4 and 8 byte edge for portable struct sizes

#define NFLOATS 5
#define NINTS 5

struct LayerStruct {
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;			// Reciprocal mu_total [cm]
	float mua;			// Absorption coefficient [1/cm]
	float g;			// Anisotropy factor [-]
	float n;			// Refractive index [-]
};

struct BoundaryStruct {
	uint32_t isHeightfield;
	uint32_t n;
	float* heights;
	float* spacings;
};

struct SimulationStruct {
	uint32_t number_of_photons;
	uint32_t ignoreAdetection;
	uint32_t n_layers;
	uint32_t start_weight;
	uint32_t begin,end; 		// mci file position offsets
	char outp_filename[STR_LEN];
	char inp_filename[STR_LEN];
	char AorB, padding[7]; 	// enforce 8 byte alignment
	struct DetStruct {
		float dr;				// Detection grid resolution, r-direction [cm]
		float dz;				// Detection grid resolution, z-direction [cm]
		uint32_t na;			// Number of grid elements in angular-direction [-]
		uint32_t nr;			// Number of grid elements in r-direction
		uint32_t nz;			// Number of grid elements in z-direction
		uint32_t padding; 		// enforce 8 byte alignment
	} det;
	LayerStruct* layers;
	struct BoundaryStruct* boundaries;
};

#endif //CUDAMCML_H