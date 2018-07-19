// uncomment to use double if supported
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// uncomment to use 64 bit atomics if supported
//#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

//#include "random.cl" // not portable since this file can end up in temp folder

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

float henyeyGreenstein(float g, float rand) {
	if (g != 0.0f) {
		return (1.0f / 2.0f * g) * (1 + g * g - pow((1 - g * g) / (1 - g + 2 * g * rand), 2));
	} else {
		return 2 * rand - 1;
	}
}

#define MAX_ITERATIONS 10000 //TODO adapt to watchdog timer

__kernel void mcml(float nAbove, float nBelow, __global struct Layer* layers, int layerCount,
volatile __global uint* reflectCount, volatile __global uint* transmitCount, volatile __global uint* absorbCount
/*float delta_r, float delta_z, float delta_a, int size_r, int size_z, int size_a,
__global __write_only float* A_rz, __global __write_only float* R_ra*/) {
	// Reflectance (specular):
	// percentage of light leaving at surface without any interaction
	// using Fesnel approximation by Schlick (no incident angle, no polarization)
	//TODO calc this on host since every thread would calc the same
	float R_specular = pow((nAbove - layers[0].n), 2) / pow((nAbove + layers[0].n), 2);
	// Q: why are the ^2 different than in Fesnel?
	// Q: is diffuse reflectance given by photons escaping at top after simulation?
	float photonWeight = 1.0f - R_specular;
	uint rng_state = wang_hash(get_global_id(0));
	float3 pos = (float3)(0.0f, 0.0f, 1.0f);
	float3 dir = (float3)(0.0f, 0.0f, 1.0f);
	int layerIndex = 0;
	for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
		__global struct Layer* currentLayer = &layers[layerIndex];
		float interactCoeff = currentLayer->absorbCoeff + currentLayer->scatterCoeff;
		// randomize step length
		rng_state = rand_xorshift(rng_state);
		float rand = (float)rng_state * (1.0f / 4294967296.0f);
		float s = -log(rand) / interactCoeff; // (noted that first s for first thread becomes infinity with current rng)
		// test layer interaction by intersection
		__global struct Boundary* intersectedBoundary = 0;
		int otherLayerIndex = 0;
		float otherN = 0;
		volatile __global uint* photonCounter = 0;
		float pathLenToIntersection = intersect(pos, dir, currentLayer->top);
		if (pathLenToIntersection >= 0 && pathLenToIntersection <= s) {
			intersectedBoundary = &currentLayer->top;
			otherLayerIndex = layerIndex - 1;
			otherN = otherLayerIndex < 0 ? nAbove : layers[otherLayerIndex].n;
			photonCounter = reflectCount;
		} else if ((pathLenToIntersection = intersect(pos, dir, currentLayer->bottom)) >= 0 && pathLenToIntersection <= s) {
			intersectedBoundary = &currentLayer->bottom;
			otherLayerIndex = layerIndex + 1;
			otherN = otherLayerIndex >= layerCount ? nBelow : layers[otherLayerIndex].n;
			photonCounter = transmitCount;
		}
		if (intersectedBoundary) {
			pos += dir * pathLenToIntersection; // unfinished part of s can be ignored
			// decide transmit or reflect
			float3 normal = (float3)(intersectedBoundary->nx, intersectedBoundary->ny, intersectedBoundary->nz);
			float cosIncident = dot(normal, -dir);
			float incidentAngle = acos(cosIncident);
			float sinTransmit = currentLayer->n * sin(incidentAngle) / otherN; // Snell's law
			float transmitAngle = asin(sinTransmit);
			float fresnelR = 1.0f/2.0f * (pow(sin(incidentAngle - transmitAngle), 2) / pow(sin(incidentAngle + transmitAngle), 2) + pow(tan(incidentAngle - transmitAngle), 2) / pow(tan(incidentAngle + transmitAngle), 2));
			rng_state = rand_xorshift(rng_state);
			float rand = (float)rng_state * (1.0f / 4294967296.0f);
			bool reflect = rand <= fresnelR;
			if (!reflect) {
				layerIndex = otherLayerIndex;
				if (layerIndex < 0 || layerIndex >= layerCount) {
					atomic_add(photonCounter, 1u);
					break;
				}
			} else {
				dir += normal * dot(normal, dir) * 2.0f; // mirror dir vector against boundary plane
			}
		} else { // absorb and scatter
			pos += dir * s; // hop
			photonWeight -= photonWeight * currentLayer->absorbCoeff / interactCoeff; // drop
			if (photonWeight < 0.0001f) {
				// roulette
				rng_state = rand_xorshift(rng_state);
				float rand = (float)rng_state * (1.0f / 4294967296.0f);
				if (rand <= 1.0f/10.0f) {
					photonWeight *= 10.0f;
				} else {
					photonWeight = 0.0f;
					atomic_add(absorbCount, 1u);
					break;
				}
			}
			// spin
			rng_state = rand_xorshift(rng_state);
			float rand = (float)rng_state * (1.0f / 4294967296.0f);
			float cosTheta = henyeyGreenstein(currentLayer->g, rand);
			float theta = acos(cosTheta);
			rng_state = rand_xorshift(rng_state);
			rand = (float)rng_state * (1.0f / 4294967296.0f);
			float psi = 2 * 3.14159265359f * rand;
			dir.x = (sin(theta) / sqrt(1.0f - dir.z * dir.z)) * (dir.x * dir.z * cos(psi) - dir.y * sin(psi)) + dir.x * cos(theta);
			dir.y = (sin(theta) / sqrt(1.0f - dir.z * dir.z)) * (dir.y * dir.z * cos(psi) - dir.x * sin(psi)) + dir.y * cos(theta);
			dir.z = -sin(theta) * cos(psi) * sqrt(1.0f - dir.z * dir.z) + dir.z * cos(theta);
			// Q: why different formula for dir close to normal?
		}
	}
}

// Basic monte carlo photon transport walkthrough
// https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/monte-carlo-methods-in-practice/monte-carlo-simulation

// Functions to get work group info
// https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/workItemFunctions.html

// OpenCL 2.0 atomics with float support etc.
// https://software.intel.com/en-us/articles/using-opencl-20-atomics#_Toc398048807