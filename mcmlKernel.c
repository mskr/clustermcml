/*********************************************************************************
*
* This code performs the MCML photon transport and is executed on the GPU.
*
*********************************************************************************/

// uncomment to use double if supported
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// uncomment to use 64 bit atomics if supported
//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

#ifdef NO_GPU
#include "randomlib.h" // rand_lcg, rand_xorshift
#else
// Include implementation as this is quicker
// than fiddling with linking in OpenCL compiler.
#include "randomlib.c"
#endif

#include "classert.h" // DEBUG_BUFFER_ARG, classert

#include "clmem.h" // CLMEM_ACCESS_ARRAY, CLMEM_ACCESS_ARRAY2D, CLMEM_ACCESS_AOS

#include "Boundary.h"
#include "Layer.h"
#include "PhotonTracker.h"

#include "geometrylib.c"

#define BoundaryArray __global struct Boundary*
#define LayerArray __global struct Layer*
#define Weight ulong
#define WeightArray volatile __global Weight*

#undef PI
#define PI 3.14159265359f

/**
* get float between 0 and 1
* just a wrapper function that picks a method from randomlib
*/
float getRandom(uint* rng_state) {

	#ifdef NO_GPU

	// use RNG from original MCML, so that our NO_GPU build guarantees exact same behavior
	return MCMLRandomNum();

	#else

	// use xorshift because it is simple and is considered to have very high quality
	// prevent being stuck with xorshift(0)==0 by fallback to lcg
	*rng_state = (*rng_state > 0) ? rand_xorshift(*rng_state) : rand_lcg(*rng_state);
	return (float)(*rng_state) * RAND_NORM;

	#endif
}

/**
* return cos of angle, which is
* more probable to be small the greater g is and
* evenly distributed if g == 0
*/
float sampleHenyeyGreenstein(uint* rng_state, float g) {
	float rand = getRandom(rng_state);
	if (g != 0.0f) {
		return (1.0f / (2.0f * g)) * (1 + g * g - pow((1 - g * g) / (1 - g + 2 * g * rand), 2));
	} else {
		return 2 * rand - 1;
	}
}

/**
* atomically add 32 bit unsigned integer to 64 bit unsigned integer
*/
void add(volatile __global ulong* dst64, uint src32) {
	// First try to add to least significant half
	// If there was an overflow, add 1 to most significant half
	if (atomic_add((volatile __global uint*)dst64, src32) + src32 < src32) {
		atomic_add(((volatile __global uint*)dst64) + 1, 1u);
	}
}

/**
* return new photon direction
* theta: angle to original direction
* psi: position on circle around original direction
*/
float3 spin(float3 dir, float theta, float psi) {
	float3 tmp = (float3)(0,0,0);
	if (fabs(dir.z) > 0.99999) {
		// when photon travels straight down the z axis
		// the regular coordinate transform would set x and y direction to (nearly) zero
		// for this case exists a simplified equivalent formula
		tmp.x = sin(theta) * cos(psi);
		tmp.y = sin(theta) * sin(psi);
		tmp.z = sign(dir.z) * cos(theta);
	} else {
		// spherical to cartesian coordinates
		tmp.x = (sin(theta) / sqrt(1.0f - dir.z * dir.z)) * (dir.x * dir.z * cos(psi) - dir.y * sin(psi)) + dir.x * cos(theta);
		tmp.y = (sin(theta) / sqrt(1.0f - dir.z * dir.z)) * (dir.y * dir.z * cos(psi) + dir.x * sin(psi)) + dir.y * cos(theta);
		tmp.z = -sin(theta) * cos(psi) * sqrt(1.0f - dir.z * dir.z) + dir.z * cos(theta);
	}
	return tmp;
}

/**
* return if photon is killed and update its weight (0 == dead)
*/
bool roulette(uint* rng_state, float* photonWeight) {
	float rand = getRandom(rng_state);
	if (rand <= 0.1f) {
		*photonWeight *= 10.0f;
	} else {
		*photonWeight = 0;
		return true;
	}
	return false;
}

/**
* Return index difference from current to adjacent layer in case
* there is an intersection with the boundary in the photon step.
* Also outputs collision normal and path length to intersection.
*/
int detectBoundaryCollision(int currentLayer, struct Line3 line,
	BoundaryArray boundaries, __global float* heights, __global float* spacings,
	float3* normal, float* pathLenToIntersection) {

	// Find intersection with top boundary
	if (boundaries[currentLayer].isHeightfield) {
		*pathLenToIntersection = intersectHeightfield(line, boundaries[currentLayer].heightfield, heights, spacings, normal);
	} else {
		const float3 middle = (float3)(0.0f, 0.0f, boundaries[currentLayer].z);
		*normal = (float3)(0.0f, 0.0f, 1.0f);
		struct Plane3 plane = {middle, *normal};
		*pathLenToIntersection = intersectPlaneWithLine(plane, line);
	}

	if (*pathLenToIntersection >= 0) {
		return -1;
	}

	// Find intersection with bottom boundary
	if (boundaries[currentLayer+1].isHeightfield) {
		*pathLenToIntersection = intersectHeightfield(line, boundaries[currentLayer+1].heightfield, heights, spacings, normal);
	} else {
		const float3 middle = (float3)(0.0f, 0.0f, boundaries[currentLayer+1].z);
		*normal = (float3)(0.0f, 0.0f, -1.0f);
		struct Plane3 plane = {middle, *normal};
		*pathLenToIntersection = intersectPlaneWithLine(plane, line);
	}

	if (*pathLenToIntersection >= 0) {
		return 1;
	}
	return 0;
}







/**
* update current layer and return if photon left the simulation domain
* the corresponding detection array is also updated
*/
bool transmit(float3 pos, float3* dir, float transmitAngle, float cosIncident, float n1, float n2, float3 normal,
float* photonWeight, int* currentLayer, int otherLayer, int layerCount,
int size_r, int size_a, float delta_r,
WeightArray R_ra, WeightArray T_ra) {

	*currentLayer = otherLayer;

	if (*currentLayer < 0) {

		// photon escaped at top => record diffuse reflectance
		// calc indices r,a
		float r = length(pos.xy);
		int r_i = (int)floor(r / delta_r);
		r_i = min(r_i, size_r - 1); // all overflowing values are accumulated at the edges

		int a_i = (int)floor(transmitAngle*2.f/PI*size_a); // using the CUDAMCML calculation here, 
		// since we are also using their IO (apparantly they convert to degrees during IO)

		add(&CLMEM_ACCESS_ARRAY2D(R_ra, __global Weight, size_r, a_i, r_i), (uint)(*photonWeight * 0xFFFFFFFF));
		
		// photon is terminated
		*photonWeight = 0;

		return true;

	} else if (*currentLayer >= layerCount) {
		
		// photon escaped at bottom => record transmittance
		// calculation analog to reflectance above
		float r = length(pos.xy);
		int r_i = (int)floor(r / delta_r);
		r_i = min(r_i, size_r - 1);
		
		int a_i = (int)floor(transmitAngle*2.f/PI*size_a);

		add(&CLMEM_ACCESS_ARRAY2D(T_ra, __global Weight, size_r, a_i, r_i), (uint)(*photonWeight * 0xFFFFFFFF));
		
		*photonWeight = 0;

		return true;
	}

	// update direction according to refraction
	if (normal.x!=(*dir).x || normal.y!=(*dir).y || normal.z!=(*dir).z) {

		float angOut = asin(n1/n2*sin(acos(cosIncident))); // Snell's law
        //FIXME: refraction was already calculated in decideReflectOrTransmit => may have wrong results when testing multiple layers!

		*dir = rotatePointAroundVector(acos(cosIncident)-angOut, cross(normal, *dir), *dir);
	}

	return false;

	//TODO compare perf against CUDAMCML, which performs buffer store only when bin changes
	// and uses local variable to add up the weights when same bin is accessed often
}







/**
* update photon direction
*/
void reflect(float3* dir, float3 normal) {

	// mirror dir vector against surface normal
	*dir = *dir - 2.0f * normal * dot(*dir, normal);
}






/**
* return if photon is reflected at boundary
*/
bool decideReflectOrTransmit(uint* rng_state, float3 dir,
LayerArray layers, int currentLayer, int otherLayer, int layerCount, float nAbove, float nBelow,
float3 normal, bool topOrBottom, float* outTransmitAngle, float* outCosIncident, float* outN1, float* outN2) {

	float otherN = otherLayer < 0 ? nAbove : otherLayer >= layerCount ? nBelow : layers[otherLayer].n;
	float cosIncident = dot(normal, -dir);

	// compute fresnel reflectance which depends on incident and transmission/refraction angles
	float fresnelR, transmitAngle;

	// straight transmission if refractive index is const
	if (otherN == layers[currentLayer].n) {
		*outTransmitAngle = 0;
		*outCosIncident = cosIncident;
		*outN1 = layers[currentLayer].n;
		*outN2 = otherN;
		return false;
	}

	// shortcut
	if (layers[currentLayer].n < otherN && otherN * otherN * (1.0f - cosIncident * cosIncident)) {
		return true;
	}

	// compute refraction
	float incidentAngle = acos(cosIncident);
	float sinTransmit = layers[currentLayer].n * sin(incidentAngle) / otherN; // Snell's law

	// more shortcuts
	float cos_crit0 = layers[currentLayer].n > otherN ? sqrt(1.0f - otherN*otherN/(layers[currentLayer].n*layers[currentLayer].n)) : 0.0f;
	float cos_crit1 = layers[currentLayer].n > otherN ? sqrt(1.0f - otherN*otherN/(layers[currentLayer].n*layers[currentLayer].n)) : 0.0;
	if (topOrBottom && cosIncident <= cos_crit0) {
		fresnelR = 1.0f;
	} else if (cosIncident <= cos_crit1) {
		fresnelR = 1.0f;
	}
	else if (cosIncident == 1.0f) {
		fresnelR = (layers[currentLayer].n - otherN) / (layers[currentLayer].n + otherN);
		fresnelR *= fresnelR;
	}
	else if (sinTransmit >= 1.0f) {
		transmitAngle = PI/2.0f;
		fresnelR = 1.0f;
	} else {

		// no valid shortcuts, calculate whole fresnel formula
		transmitAngle = asin(sinTransmit);
		fresnelR = 1.0f/2.0f * (pow(sin(incidentAngle - transmitAngle), 2) / pow(sin(incidentAngle + transmitAngle), 2) + pow(tan(incidentAngle - transmitAngle), 2) / pow(tan(incidentAngle + transmitAngle), 2));
	}

	// output info to update photon direction
	*outTransmitAngle = transmitAngle;
	*outCosIncident = cosIncident;
	*outN1 = layers[currentLayer].n;
	*outN2 = otherN;

	// photon is reflected if random threshold exceeded (according to mcml)
	float rand = getRandom(rng_state);
	return (rand <= fresnelR);
}




// MAX_ITERATIONS of main loop, i.e. photon bounces, to control time
// spent on the GPU per round.
// After each round, the CPU checks all photons for their weight and
// restarts a photon if it finished and total target number not reached.
// If MAX_ITERATIONS is...
// ... small, more GPU rounds and more weight-checks are needed, but
// also finished photons can be restarted earlier, which reduces
// inactive threads.
// ... large, process might be killed by Windows, or the first call to
// clEnqueueReadBuffer might fail with CL_OUT_OF_RESOURCES.
#define MAX_ITERATIONS 1000





/**
* Monte carlo photon transport in multilayered media.
* 
* Input are arrays of layers and boundaries as well as
* the diffractive indices of the media above and below.
* 
* Output is consisting of
* Reflectance (diffuse), Absorption and Transmittance,
* detected on a grid whose one dimension is
* radius, meaning distance from origin in xy-plane
* and other dimension is escape angle resp. z-depth.
* The grid has configurable size and bin delta.
* Delta units are centimeters.
* Angle delta is implicitly given as angles will be
* between 0 and 90 degrees.
*
* This is a kernel function able to run on GPUs.
* Using the GPU for too long would cause the OS
* to terminate the process. Therefore a fixed
* number of photon bounces is done per round.
* To keep track of photon states across rounds,
* a tracker array is maintained which size equals
* the number of photons that the GPU can simulate
* in parallel, i.e. the thread count.
*
* Model of boundaries and layers:
*
* --------------------------------- boundary 0
*              layer 0
* --------------------------------- boundary 1
*               ...
*/
__kernel void mcml(
float nAbove, // refractive index above
float nBelow, // refractive index below
LayerArray layers, // layer array
int layerCount, // size of layer array
BoundaryArray boundaries, // boundary array
__global float* heights, // height samples of all boundaries
__global float* spacings, // sample spacings of all boundaries
int size_r, // number of radial bins
int size_a, // number of angular bins
int size_z, // number of depth bins
float delta_r, // spacing between radial bins
float delta_z, // spacing between depth bins
WeightArray R_ra, // reflectance array
WeightArray T_ra, // transmittance array
WeightArray A_rz, // absorption array
__global struct PhotonTracker* photonStates // photon tracking buffer
DEBUG_BUFFER_ARG) // optional debug buffer
{
	// Get current photon state
	__global struct PhotonTracker* state = &photonStates[get_global_id(0)];

	if (state->isDead) {
		// This photon was not restarted because enough are in the pipeline
		return; // I have no work :(
	}
	uint rng_state = state->rngState;
	float photonWeight = state->weight;
	float3 pos = (float3)(state->x, state->y, state->z);
	float3 dir = (float3)(state->dx, state->dy, state->dz);
	int currentLayer = state->layerIndex;

	int disabledBoundary = -1;

	// Simulate MAX_ITERATIONS bounces
	for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
		float interactCoeff = layers[currentLayer].absorbCoeff + layers[currentLayer].scatterCoeff;

		// Randomize step length
		float rand = getRandom(&rng_state);
		float s = -log(rand) / interactCoeff;

		// Uncomment to output some lengths of the first step of a photon (DEBUG define required)
		// if (get_global_id(0) < 2048/4) // limited by debug buffer size
		// 	((__global float*)DEBUG_BUFFER)[get_global_id(0)] = s;
		// break;

		// Test intersection with top and bottom boundaries of current layer
		float3 normal = (float3)(0); float pathLenToIntersection = 0;
		struct Line3 line = { pos, pos + s*dir };
		int layerChange = detectBoundaryCollision(currentLayer, line, boundaries, heights, spacings, &normal, &pathLenToIntersection);
		if (layerChange != 0 && disabledBoundary != max(currentLayer, currentLayer+layerChange)) {

			// In theory, the photon should now be moved exactly onto the boundary. Another collision with the same boundary should thus not be possible in the next step.
			// Due to floating point inaccuracies, the photon may be moved slightly beyond or before the boundary. Therefore we need to memorize for the next step that this boundary is temporarly disabled for collision.
			disabledBoundary = max(currentLayer, currentLayer + layerChange);

			// Hop (unfinished part of s can be ignored)
			pos += dir * pathLenToIntersection;

			// MCML does not drop weight here for some reason

			// Transmit or reflect at boundary
			float transmitAngle = 0; float cosIncident = 0; float n1 = 0; float n2 = 0;
			if (decideReflectOrTransmit(&rng_state, dir, layers, currentLayer, currentLayer+layerChange, layerCount, nAbove, nBelow, normal, layerChange<0, &transmitAngle, &cosIncident, &n1, &n2)) {
				reflect(&dir, normal);
			} else {
				if (transmit(pos, &dir, transmitAngle, cosIncident, n1, n2, normal, &photonWeight, &currentLayer, currentLayer+layerChange, layerCount, size_r, size_a, delta_r, R_ra, T_ra)) {
					break;
				}
			}
		} else {

			disabledBoundary = -1;

			// Hop
			pos += dir * s;

			// Drop
			float dW = photonWeight * layers[currentLayer].absorbCoeff / interactCoeff;
			photonWeight -= dW;

			// Record A
			#ifndef IGNORE_A
			{
				float z = pos.z;
				int z_i = (int)floor(z / delta_z);
				z_i = max(0, z_i);
				z_i = min(z_i, size_z - 1);
				float r = length(pos.xy);
				int r_i = (int)floor(r / delta_r);
				r_i = min(r_i, size_r - 1); // all overflowing values are accumulated at the edges
				add(&CLMEM_ACCESS_ARRAY2D(A_rz, __global Weight, size_r, z_i, r_i), (uint)(dW * 0xFFFFFFFF));
			}
			#endif

			// Spin
			// for g==0 cosTheta is evenly distributed...
			float cosTheta = sampleHenyeyGreenstein(&rng_state, layers[currentLayer].g);
			// ... and theta has most values at PI/2, which is unintuitive but correct
			float theta = acos(cosTheta);
			rand = getRandom(&rng_state);
			float psi = 2 * PI * rand;
			dir = spin(dir, theta, psi);
			dir = normalize(dir); // normalize necessary wrt precision problems of float

			// Decide if photon dies
			if (photonWeight < 0.0001f) {
				if (roulette(&rng_state, &photonWeight)) {
					break;
				}
			}
		}
	}

	// Uncomment to count finished photons (DEBUG define required)
	// Something like this can be used to restart photons on the fly,
	// without exceeding target count, to avoid inactive threads.
	// if (photonWeight == 0) {
	// 	atomic_add((volatile __global uint*)DEBUG_BUFFER, 1); 
	// }

	// Uncomment to output spacings of heightfield boundaries
	// if (get_global_id(0) < 60)
	// 	((__global float*)DEBUG_BUFFER)[get_global_id(0)] = spacings[get_global_id(0)];

	// Uncomment to output parameters of first two boundaries
	// if (get_global_id(0) < 2) 
	// 	((__global uint*)DEBUG_BUFFER)[get_global_id(0)] = boundaries[get_global_id(0)].isHeightfield;
	// if (get_global_id(0) > 1 && get_global_id(0) < 4)
	// 	((__global float*)DEBUG_BUFFER)[get_global_id(0)] = boundaries[get_global_id(0)-2].heightfield.center.z;
	// if (get_global_id(0) > 3 && get_global_id(0) < 6)
	// 	((__global uint*)DEBUG_BUFFER)[get_global_id(0)] = boundaries[get_global_id(0)-4].heightfield.i_heights;
	// if (get_global_id(0) > 5 && get_global_id(0) < 8)
	// 	((__global uint*)DEBUG_BUFFER)[get_global_id(0)] = boundaries[get_global_id(0)-6].heightfield.n_heights;
	// if (get_global_id(0) > 7 && get_global_id(0) < 10)
	// 	((__global uint*)DEBUG_BUFFER)[get_global_id(0)] = boundaries[get_global_id(0)-8].heightfield.i_spacings;
	// if (get_global_id(0) > 9 && get_global_id(0) < 12)
	// 	((__global uint*)DEBUG_BUFFER)[get_global_id(0)] = boundaries[get_global_id(0)-10].heightfield.n_spacings;
	// if (get_global_id(0) > 11 && get_global_id(0) < 13)
	// 	((__global uint*)DEBUG_BUFFER)[get_global_id(0)] = sizeof(boundaries[get_global_id(0)]);

	// Save state for the next round
	state->x = pos.x; state->y = pos.y; state->z = pos.z;
	state->dx = dir.x; state->dy = dir.y; state->dz = dir.z;
	state->weight = photonWeight;
	state->layerIndex = currentLayer;
	state->rngState = rng_state;
}

// Basic monte carlo photon transport walkthrough
// https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/monte-carlo-methods-in-practice/monte-carlo-simulation

// Functions to get work group info
// https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/workItemFunctions.html

// OpenCL atomics
// https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/atomicFunctions.html

// OpenCL 2.0 atomics with float support etc.
// https://software.intel.com/en-us/articles/using-opencl-20-atomics#_Toc398048807