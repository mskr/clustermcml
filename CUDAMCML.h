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

#define STR_LEN 200

#define NFLOATS 5
#define NINTS 5

// TYPEDEFS
typedef struct 
{
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;			// Reciprocal mu_total [cm]
	float mua;			// Absorption coefficient [1/cm]
	float g;			// Anisotropy factor [-]
	float n;			// Refractive index [-]
}LayerStruct;

typedef struct 
{
	float dr;		// Detection grid resolution, r-direction [cm]
	float dz;		// Detection grid resolution, z-direction [cm]
	
	int na;			// Number of grid elements in angular-direction [-]
	int nr;			// Number of grid elements in r-direction
	int nz;			// Number of grid elements in z-direction

}DetStruct;


typedef struct 
{
	unsigned long number_of_photons;
	int ignoreAdetection;
	unsigned int n_layers;
	unsigned int start_weight;
	char outp_filename[STR_LEN];
	char inp_filename[STR_LEN];
	long begin,end;
	char AorB;
	DetStruct det;
	LayerStruct* layers;
}SimulationStruct;

#endif //CUDAMCML_H