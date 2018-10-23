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

#undef PI
#define PI 3.14159265359f

/**
* return path length to ray-plane intersection
*/
float intersectPlane(float3 pos, float3 dir, float3 point, float3 normal) {
	float a = dot((point - pos), normal);
	if (a > -1e-6f) return INFINITY; // behind plane
	float b = dot(dir, normal);
	if (b > -1e-6f) return INFINITY; // facing away
	float pathLenToIntersection = a / b;
	return pathLenToIntersection;
}

/**
* Return cartesian z value from radial heightfield at cartesian position xy.
* Also output normal that is oriented in -z direction.
* Heights are interpreted in -z direction since +z goes down into material.
*/
float readRadialHeightfield(float2 pos, float3 center, float heightfield[BOUNDARY_SAMPLES], float3* outNormal) {
	float res = (float)BOUNDARY_WIDTH / (float)BOUNDARY_SAMPLES;
	// Calc p in heightmap coordinate system
	float2 p = pos.xy - center.xy;
	// Calc radial offset
	float r = length(p);
	// Calc array offset using sampling resolution
	float x = r / res;
	int i = (int)x;
	// Get 2 nearest samples
	float h0 = heightfield[min(i, BOUNDARY_SAMPLES-1)];
	float h1 = heightfield[min(i+1, BOUNDARY_SAMPLES-1)];
	// Calc distances from the 2 samples
	float f = fract(x); float rf = 1.0f-f;
	// Calc gradient from p1 to p0,
	// so that p0 lies on the nearest sample towards center
	// and p1 lies on the next outer sample
	float2 p0 = (float2)(0, 0);
	float2 p1 = (float2)(res, 0);
	if (x > 0) {
		p0 = p*(1.0f-f/x);
		p1 = p*(1.0f+rf/x);
	}
	float3 gradient = normalize((float3)(p0, -h0) - (float3)(p1, -h1));
	// Calc clockwise tangent
	float3 tangent = normalize((float3)(p.y, -p.x, 0));
	// Calc normal from gradient and tangent
	// note: cross(a,b) forms right-handed system, where a==thumb
	// then translate normal back into cartesian coordinates
	*outNormal = normalize(cross(tangent, gradient) + center);
	// Linear interpolation
	return center.z - mix(h0, h1, f);
}

/**
* find intersection of line with radial heightfield
*/
float intersectHeightfield(float3 pos, float3 dir, float len, float3 center, float heightfield[BOUNDARY_SAMPLES], float3* outNormal) {
	// Coarse raymarching search
	float D = 0.00001f;
	float3 p = pos;
	float lastDz = p.z - readRadialHeightfield(p.xy, center, heightfield, outNormal);
	float d = D;
	float3 ds = D * dir;
	while (d < len) {
		p += ds;
		float dz = p.z - readRadialHeightfield(p.xy, center, heightfield, outNormal);
		if (sign(dz) != sign(lastDz)) {
			// if ray came from +z, i.e. (dz < 0), normal will point down
			if (dz < 0) *outNormal *= -1.0f;
			// Path to intersection
			return d;
		}
		lastDz = dz;
		d += D;
	}
	return -1.f;
}

//TODO try exact heightfield intersection using intermediate planes, like Maisch pointed out

/**
* return cos of angle, which is
* more probable to be small the greater g is and
* evenly distributed if g == 0
*/
float sampleHenyeyGreenstein(uint* rng_state, float g) {
	float rand = (float)(*rng_state = rand_xorshift(*rng_state)) * RAND_NORM;
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
	float rand = (float)(*rng_state = rand_xorshift(*rng_state)) * RAND_NORM;
	if (rand <= 0.1f) {
		*photonWeight *= 10.0f;
	} else {
		*photonWeight = 0;
		return true;
	}
	return false;
}

/**
* return if there is an intersection in the photon step and write according parameters
*/
bool detectBoundaryCollision(float3 pos, float3 dir, float s, int currentLayer, __global struct Boundary* boundaries, // in
__global struct Boundary** intersectedBoundary, float3 normal, int* otherLayer, float* pathLenToIntersection, bool* topOrBottom) { // out
	//TODO return -1, 0 or 1 to indicate layer change and get rid of otherLayer and topOrBottom args
	// Find intersection with top boundary
	float3 n = (float3)(0);
	*pathLenToIntersection = intersectHeightfield(pos, dir, s, (float3)(0,0,boundaries[currentLayer]), boundaries[currentLayer].heightfield, &n);
	if (*pathLenToIntersection >= 0) {
		*intersectedBoundary = &boundaries[currentLayer];
		*normal = n;
		*otherLayer = currentLayer - 1;
		*topOrBottom = true;
		return true;
	}
	// Find intersection with bottom boundary
	*pathLenToIntersection = intersectHeightfield(pos, dir, s,(float3)(0,0,boundaries[currentLayer+1]), boundaries[currentLayer+1].heightfield &n);
	if (*pathLenToIntersection >= 0) {
		*intersectedBoundary = &boundaries[currentLayer+1];
		*normal = n;
		*otherLayer = currentLayer + 1;
		*topOrBottom = false;
		return true;
	}
	return false;
}

/**
* update current layer and return if photon left the simulation domain
* the corresponding detection array is also updated
*/
bool transmit(float3 pos, float3* dir, float transmitAngle, float cosIncident, float n1, float n2,
float* photonWeight, int* currentLayer, int otherLayer, int layerCount,
int size_r, int size_a, float delta_r,
volatile __global ulong* R_ra, volatile __global ulong* T_ra) {
	*currentLayer = otherLayer;
	if (*currentLayer < 0) {
		// photon escaped at top => record diffuse reflectance
		// calc indices r,a
		float r = length(pos.xy);
		int i = (int)floor(r / delta_r);
		i = min(i, size_r - 1); // all overflowing values are accumulated at the edges
		float a = transmitAngle / (2.0f * PI) * 360.0f;
		int j = (int)floor(a / (90.0f / size_a));
		add(&R_ra[i * size_a + j], (uint)(*photonWeight * 0xFFFFFFFF));
		// photon is terminated
		*photonWeight = 0;
		return true;
	} else if (*currentLayer >= layerCount) {
		// photon escaped at bottom => record transmittance
		// calc indices r,a
		float r = length(pos.xy);
		int i = (int)floor(r / delta_r);
		i = min(i, size_r - 1); // all overflowing values are accumulated at the edges
		float a = transmitAngle / (2.0f * PI) * 360.0f;
		int j = (int)floor(a / (90.0f / size_a));
		add(&T_ra[i * size_a + j], (uint)(*photonWeight * 0xFFFFFFFF));
		*photonWeight = 0;
		return true;
	}
	// update direction
	float r = n1 / n2;
	float e = r * r * (1.0f - cosIncident * cosIncident);
	(*dir).x = r;
	(*dir).y = r;
	//TODO check if this works for all normals!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	(*dir).z = copysign(sqrt(1.0f - e), (*dir).z);
	return false;

	//TODO compare perf against CUDAMCML, which performs buffer store only when bin changes
	// and uses local variable to add up the weights when same bin is accessed often
}

/**
* update photon direction
*/
void reflect(float3* dir, __global struct Boundary* intersectedBoundary) {
	float3 normal = (float3)(intersectedBoundary->nx, intersectedBoundary->ny, intersectedBoundary->nz);
	// mirror dir vector against boundary plane
	*dir = *dir - 2.0f * normal * dot(*dir, normal);
}

/**
* return if photon is reflected at boundary
*/
bool decideReflectOrTransmit(uint* rng_state, float3 dir,
__global struct Layer* layers, int currentLayer, int otherLayer, int layerCount, float nAbove, float nBelow,
__global struct Boundary* intersectedBoundary, bool topOrBottom, float* outTransmitAngle, float* outCosIncident, float* outN1, float* outN2) {
	float otherN = otherLayer < 0 ? nAbove : otherLayer >= layerCount ? nBelow : layers[otherLayer].n;
	float3 normal = (float3)(intersectedBoundary->nx, intersectedBoundary->ny, intersectedBoundary->nz);
	float cosIncident = dot(normal, -dir);
	float fresnelR, transmitAngle;
	// straight transmission if refractive index is const
	if (otherN == layers[currentLayer].n) {
		*outTransmitAngle = 0;
		*outCosIncident = cosIncident;
		*outN1 = layers[currentLayer].n;
		*outN2 = otherN;
		return false;
	}
	if (layers[currentLayer].n < otherN && otherN * otherN * (1.0f - cosIncident * cosIncident)) {
		return true;
	}
	float incidentAngle = acos(cosIncident);
	float sinTransmit = layers[currentLayer].n * sin(incidentAngle) / otherN; // Snell's law
	float cos_crit0 = layers[currentLayer].n > otherN ? sqrt(1.0f - otherN*otherN/(layers[currentLayer].n*layers[currentLayer].n)) : 0.0f;
	float cos_crit1 = layers[currentLayer].n > otherN ? sqrt(1.0 - otherN*otherN/(layers[currentLayer].n*layers[currentLayer].n)) : 0.0;
	if (topOrBottom && cosIncident <= cos_crit0) {
		fresnelR = 1.0f;
	} else if (cosIncident <= cos_crit1) {
		fresnelR = 1.0f;
	}
	//TODO We still transmit too much photons! Find a definite state where mcml does not transmit but we do (just output the index of the photon).
	else if (cosIncident == 1.0f) {
		fresnelR = (layers[currentLayer].n - otherN) / (layers[currentLayer].n + otherN);
		fresnelR *= fresnelR;
	}
	else if (sinTransmit >= 1.0f) {
		transmitAngle = PI/2.0f;
		fresnelR = 1.0f;
	} else {
		transmitAngle = asin(sinTransmit);
		fresnelR = 1.0f/2.0f * (pow(sin(incidentAngle - transmitAngle), 2) / pow(sin(incidentAngle + transmitAngle), 2) + pow(tan(incidentAngle - transmitAngle), 2) / pow(tan(incidentAngle + transmitAngle), 2));
	}
	float rand = (float)(*rng_state = rand_xorshift(*rng_state)) * RAND_NORM;
	*outTransmitAngle = transmitAngle;
	*outCosIncident = cosIncident;
	*outN1 = layers[currentLayer].n;
	*outN2 = otherN;
	return (rand <= fresnelR);
}




// control time spent on the GPU in each round
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
*/
__kernel void mcml(
float nAbove, // refractive index above
float nBelow, // refractive index below
__global struct Layer* layers, // layer array
int layerCount, // size of layer array
__global struct Boundary* boundaries, // boundary array
int size_r, // number of radial bins
int size_a, // number of angular bins
int size_z, // number of depth bins
float delta_r, // spacing between radial bins
float delta_z, // spacing between depth bins
volatile __global ulong* R_ra, // reflectance array
volatile __global ulong* T_ra, // transmittance array
volatile __global ulong* A_rz, // absorption array
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
	if (rng_state == 0) rng_state = wang_hash(get_global_id(0));
	float photonWeight = state->weight;
	float3 pos = (float3)(state->x, state->y, state->z);
	float3 dir = (float3)(state->dx, state->dy, state->dz);
	int currentLayer = state->layerIndex;

	// Simulate a few bounces
	for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
		float interactCoeff = layers[currentLayer].absorbCoeff + layers[currentLayer].scatterCoeff;

		// Randomize step length
		// prevent being stuck with xorshift(0)==0 by fallback to lcg
		rng_state = (rng_state > 0) ? rand_xorshift(rng_state) : rand_lcg(rng_state);
		float rand = (float)rng_state * RAND_NORM;
		float s = -log(rand) / interactCoeff;

		// Uncomment to output some lengths of the first step of a photon (DEBUG mode required)
		// if (get_global_id(0) < 2048/4) // limited by debug buffer size
		// 	((__global float*)DEBUG_BUFFER)[get_global_id(0)] = s;
		// break;

		// Test intersection with top and bottom boundaries of current layer
		__global struct Boundary* intersectedBoundary = 0; float3 normal = (float3)(0); int otherLayer = 0; float pathLenToIntersection = 0; bool topOrBottom = 0;
		if (detectBoundaryCollision(pos, dir, s, currentLayer, boundaries, &intersectedBoundary, &normal, &otherLayer, &pathLenToIntersection, &topOrBottom)) {

			// Hop (unfinished part of s can be ignored)
			pos += dir * pathLenToIntersection;

			// MCML does not drop weight here for some reason

			// Transmit or reflect at boundary
			float transmitAngle = 0; float cosIncident = 0; float n1 = 0; float n2 = 0;
			if (decideReflectOrTransmit(&rng_state, dir, layers, currentLayer, otherLayer, layerCount, nAbove, nBelow, 
			intersectedBoundary, topOrBottom, &transmitAngle, &cosIncident, &n1, &n2)) {
				reflect(&dir, intersectedBoundary);
			} else {
				if (transmit(pos, &dir, transmitAngle, cosIncident, n1, n2, &photonWeight, &currentLayer, otherLayer, layerCount,
				size_r, size_a, delta_r, R_ra, T_ra)) {
					break;
				}
			}
		} else {

			// Hop
			pos += dir * s;

			// Drop
			float dW = photonWeight * layers[currentLayer].absorbCoeff / interactCoeff;
			photonWeight -= dW;

			// Record A
			#ifndef IGNORE_A
			{
				float z = pos.z;
				int i = (int)floor(z / delta_z);
				i = min(i, size_z - 1);
				float r = length(pos.xy);
				int j = (int)floor(r / delta_r);
				j = min(j, size_r - 1); // all overflowing values are accumulated at the edges
				add(&A_rz[i * size_r + j], (uint)(dW * 0xFFFFFFFF));
			}
			#endif

			// Spin
			// for g==0 cosTheta is evenly distributed...
			float cosTheta = sampleHenyeyGreenstein(&rng_state, layers[currentLayer].g);
			// ... and theta has most values at PI/2, which is unintuitive but correct
			float theta = acos(cosTheta);
			rand = (float)(rng_state = rand_xorshift(rng_state)) * RAND_NORM;
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