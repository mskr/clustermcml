// uncomment to use double if supported
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// uncomment to use 64 bit atomics if supported
//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

//#include "random.cl" // not portable since this file can end up in temp folder

#define PI 3.14159265359f

// An assert macro that writes error message to host buffer and returns from current function
//TODO also write values of local variables or dump the whole stack frame
#ifdef DEBUG
#define DEBUG_BUFFER_ARG ,__global char* debugBuffer
#define STR_COPY(src, dst) for(int i=0; src[i]!='\0';i++) dst[i]=src[i];
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x) //The extra level of indirection will allow the preprocessor to expand the macros before they are converted to strings.
#define assert(expr, t, a, b, c)\
	if(!(expr)) {\
		const __constant char* msg = "error: assertion failed at line " STR(__LINE__);\
		int i=0; for(; msg[i]!='\0';i++) debugBuffer[i]=msg[i];\
		((__global t*)debugBuffer)[i] = a; ((__global t*)debugBuffer)[i+1] = b; ((__global t*)debugBuffer)[i+2] = c;\
		return;\
	}
#else
#define DEBUG_BUFFER_ARG
#define assert(expr, t, a, b, c)
#endif

// Following Nathan Reed's article:
// http://reedbeta.com/blog/quick-and-easy-gpu-random-numbers-in-d3d11/
// PRNGs are designed to go deep, i.e. have good distributions when sequentially updating state
// Hashes are designed to go wide, i.e. have good distributions across initial seeds
// Using thread index as seed, hashes map better to the GPU
// Hashed thread index can also be used as seed for the PRNGs

// For normalized random number in [0, 1) use: 
// (float)rng_state * (1.0f / 4294967296.0f)

// Xorshift algorithm from George Marsaglia's paper
uint rand_xorshift(uint rng_state) {
	rng_state ^= (rng_state << 13);
	rng_state ^= (rng_state >> 17);
	rng_state ^= (rng_state << 5);
	return rng_state;
}

// LCG (Linear congruential generator) values from Numerical Recipes
uint rand_lcg(uint rng_state) {
	rng_state = 1664525 * rng_state + 1013904223;
	return rng_state;
}

// Hash function by Thomas Wang
// http://www.burtleburtle.net/bob/hash/integer.html
uint wang_hash(uint seed) {
	seed = (seed ^ 61) ^ (seed >> 16);
	seed *= 9;
	seed = seed ^ (seed >> 4);
	seed *= 0x27d4eb2d;
	seed = seed ^ (seed >> 15);
	return seed;
}

// Mersenne Twister
// https://us.fixstars.com/opencl/book/OpenCLProgrammingBook/mersenne-twister/

// Multiply With Carry
// https://www.ast.cam.ac.uk/~stg20/cuda/random/index.html
// used by CUDAMCML, they also have versions for [0,1) and (0,1]

// Function to integrate by simpson kernel
float simpson_f(float x) {
	return 4 / (1 + x*x);
}

// Simpson integration
__kernel void simpson(float ha, float hb, __global float* out) {
	size_t i = get_global_id(0);
	size_t nthreads = get_global_size(0);
	float a = ha + (hb - ha) * i/nthreads;
	float b = ha + (hb - ha) * (i+1)/nthreads;
	float value = simpson_f(a)/2 + simpson_f(b)/2;
	float xleft;
	float x = a;
	int n = 1;
	for (int i = 1; i < n; ++i) {
		xleft = x;
		x = a + i * (b - a) / n;
		value += simpson_f(x) + 2 * simpson_f((xleft + x)/2);
	}
	value += 2 * simpson_f((x + b)/2);
	value *= (b - a) / n / 3;
	out[i] = value;
}

// Monte carlo approximation of PI
__kernel void mcpi(const int npoints, __global uint* out) {
	size_t i = get_global_id(0);
	uint rng_state = wang_hash(i);
	uint count = 0;
	for (int j = 0; j < npoints; j++) {
		rng_state = rand_xorshift(rng_state);
		float x = (float)rng_state * (1.0f / 4294967296.0f);
		rng_state = rand_xorshift(rng_state);
		float y = (float)rng_state * (1.0f / 4294967296.0f);
		// check if inside quarter unit circle
		if (x * x + y * y < 1.0f) {
			count++;
		}
	}
	out[i] = count;
}

/***** MCML *****/

// Boundary with customizable shapes
// Q: Which shapes make sense with respect to layers intersecting each other?
struct Boundary {
	float z; // depth
	float nx, ny, nz; // normal
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

float intersect(float3 pos, float3 dir, struct Boundary bound) {
	// find ray-plane intersection point
	// if found return path length to intersection
	// otherwise return negative value
	float3 normal = (float3)(bound.nx, bound.ny, bound.nz);
	float a = dot(((float3)(0.0f, 0.0f, bound.z) - pos), normal);
	if (a > -1e-6f) return -1.0f; // behind plane
	float b = dot(dir, normal);
	if (b > -1e-6f) return -1.0f; // facing away
	float pathLenToIntersection = a / b;
	return pathLenToIntersection;
}

float sampleHenyeyGreenstein(float g, float rand) {
	if (g != 0.0f) {
		return (1.0f / 2.0f * g) * (1 + g * g - pow((1 - g * g) / (1 - g + 2 * g * rand), 2));
	} else {
		return 2 * rand - 1;
	}
}

void add(volatile __global ulong* dst64, uint src32) {
	// Add 32 bit unsigned integer to 64 bit unsigned integer, atomically
	// First try to add to least significant half
	// If there was an overflow, add 1 to most significant half
	if (atomic_add((volatile __global uint*)dst64, src32) + src32 < src32) {
		atomic_add(((volatile __global uint*)dst64) + 1, 1u);
	}
}

float3 spin(float3 dir, float theta, float psi) {
	if (fabs(dir.z) > 0.99999) {
		// when photon travels straight down the z axis
		// the regular coordinate transform would set x and y direction to (nearly) zero
		dir.x = sin(theta) * cos(psi);
		dir.y = sin(theta) * sin(psi);
		dir.z = sign(dir.z) * cos(theta);
	} else {
		// spherical to cartesian coordinates
		dir.x = (sin(theta) / sqrt(1.0f - dir.z * dir.z)) * (dir.x * dir.z * cos(psi) - dir.y * sin(psi)) + dir.x * cos(theta);
		dir.y = (sin(theta) / sqrt(1.0f - dir.z * dir.z)) * (dir.y * dir.z * cos(psi) - dir.x * sin(psi)) + dir.y * cos(theta);
		dir.z = -sin(theta) * cos(psi) * sqrt(1.0f - dir.z * dir.z) + dir.z * cos(theta);
	}
	return dir;
}

bool roulette(uint* rng_state, float* photonWeight) {
	*rng_state = rand_xorshift(*rng_state);
	float rand = (float)(*rng_state) * (1.0f / 4294967296.0f);
	if (rand <= 1.0f/10.0f) {
		*photonWeight *= 10.0f;
	} else {
		*photonWeight = 0;
		return true;
	}
	return false;
}

bool findIntersection(float3 pos, float3 dir, float s, __global struct Layer* layers, int currentLayer,
__global struct Boundary** intersectedBoundary, int* otherLayer, float* pathLenToIntersection) {
	if ((*pathLenToIntersection = intersect(pos, dir, layers[currentLayer].top)) >= 0 && *pathLenToIntersection <= s) {
		*intersectedBoundary = &layers[currentLayer].top;
		*otherLayer = currentLayer - 1;
		return true;
	}
	if ((*pathLenToIntersection = intersect(pos, dir, layers[currentLayer].bottom)) >= 0 && *pathLenToIntersection <= s) {
		*intersectedBoundary = &layers[currentLayer].bottom;
		*otherLayer = currentLayer + 1;
		return true;
	}
	return false;
}

bool transmit(float3 pos, float transmitAngle, float* photonWeight, int* currentLayer, int otherLayer, int layerCount,
int size_r, int size_a, float delta_r, volatile __global ulong* R_ra) {
	//TODO update direction
	*currentLayer = otherLayer;
	if (*currentLayer < 0) { // photon escaped at top => record diffuse reflectance
		float r = length(pos.xy);
		int i = (int)floor(r / delta_r);
		i = min(i, size_r - 1); // all overflowing values are accumulated at the edges
		float a = transmitAngle / (2.0f * PI) * 360.0f;
		int j = (int)floor(a / (90.0f / size_a));
		add(&R_ra[i * size_a + j], (uint)(*photonWeight * 0xFFFFFFFF));
		*photonWeight = 0;
		return true;
	} else if (*currentLayer >= layerCount) {
		*photonWeight = 0;
		return true;
	}
	return false;
}

void reflect(float3* dir, float3 normal) {
	*dir += normal * dot(normal, *dir) * 2.0f; // mirror dir vector against boundary plane
}

#define MAX_ITERATIONS 1000

__kernel void mcml(float nAbove, float nBelow, __global struct Layer* layers, int layerCount,
int size_r, int size_a, float delta_r, 
volatile __global ulong* R_ra,
__global struct PhotonState* photonStates
DEBUG_BUFFER_ARG)
{
	__global struct PhotonState* state = &photonStates[get_global_id(0)];
	uint rng_state = wang_hash(get_global_id(0));
	float photonWeight = state->weight;
	float3 pos = (float3)(state->x, state->y, state->z);
	float3 dir = (float3)(state->dx, state->dy, state->dz);
	int currentLayer = state->layerIndex;
	for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
		float interactCoeff = layers[currentLayer].absorbCoeff + layers[currentLayer].scatterCoeff;
		// randomize step length
		rng_state = rand_xorshift(rng_state);
		float rand = (float)rng_state * (1.0f / 4294967296.0f); // rng_state can be max 0xFFFFFFFF=4294967295 => rand in [0,1)
		float s = -log(rand) / interactCoeff; // (noted that first s for first thread becomes infinity with current rng)
		__global struct Boundary* intersectedBoundary = 0; int otherLayer = 0; float pathLenToIntersection = 0;
		if (findIntersection(pos, dir, s, layers, currentLayer, &intersectedBoundary, &otherLayer, &pathLenToIntersection)) {
			pos += dir * pathLenToIntersection; // unfinished part of s can be ignored
			//TODO drop some weight here?
			// decide transmit or reflect
			float otherN = otherLayer < 0 ? nAbove : otherLayer >= layerCount ? nBelow : layers[otherLayer].n;
			float3 normal = (float3)(intersectedBoundary->nx, intersectedBoundary->ny, intersectedBoundary->nz);
			assert(fabs(length(normal) - 1.0f) < 0.001f, float, normal.x, normal.y, normal.z);
			float cosIncident = dot(normal, -dir);
			float incidentAngle = acos(cosIncident);
			float sinTransmit = layers[currentLayer].n * sin(incidentAngle) / otherN; // Snell's law
			float transmitAngle = asin(sinTransmit);
			float fresnelR = 1.0f/2.0f * (pow(sin(incidentAngle - transmitAngle), 2) / pow(sin(incidentAngle + transmitAngle), 2) + pow(tan(incidentAngle - transmitAngle), 2) / pow(tan(incidentAngle + transmitAngle), 2));
			rng_state = rand_xorshift(rng_state);
			float rand = (float)rng_state * (1.0f / 4294967296.0f);
			if (rand <= fresnelR) {
				reflect(&dir, normal);
			} else {
				if (transmit(pos, transmitAngle, &photonWeight, &currentLayer, otherLayer, layerCount, size_r, size_a, delta_r, R_ra)) break;
			}
		} else { // absorb and scatter
			pos += dir * s; // hop
			photonWeight -= photonWeight * layers[currentLayer].absorbCoeff / interactCoeff; // drop
			if (photonWeight < 0.0001f) {
				if (roulette(&rng_state, &photonWeight)) break;
			}
			// spin
			rng_state = rand_xorshift(rng_state);
			float rand = (float)rng_state * (1.0f / 4294967296.0f);
			float cosTheta = sampleHenyeyGreenstein(layers[currentLayer].g, rand); // for g==0 cosTheta is evenly distributed
			float theta = acos(cosTheta); // for g==0 theta has most values at pi/2, which is correct
			rng_state = rand_xorshift(rng_state);
			rand = (float)rng_state * (1.0f / 4294967296.0f);
			float psi = 2 * PI * rand;
			dir = spin(dir, theta, psi);
			dir = normalize(dir); // necessary wrt precision problems of float
		}
	}
	state->x = pos.x; state->y = pos.y; state->z = pos.z;
	state->dx = dir.x; state->dy = dir.y; state->dz = dir.z;
	state->weight = photonWeight;
	state->layerIndex = currentLayer;
}

// Basic monte carlo photon transport walkthrough
// https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/monte-carlo-methods-in-practice/monte-carlo-simulation

// Functions to get work group info
// https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/workItemFunctions.html

// OpenCL atomics
// https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/atomicFunctions.html

// OpenCL 2.0 atomics with float support etc.
// https://software.intel.com/en-us/articles/using-opencl-20-atomics#_Toc398048807