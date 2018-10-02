#define const __constant
#define uint32_t uint
#include "randomlib.c"
#undef const
#undef uint32_t

// Monte carlo approximation of PI
__kernel void mcpi(const int npoints, __global uint* out) {
	size_t i = get_global_id(0);
	uint rng_state = wang_hash(i);
	if (rng_state == 0) rng_state = 1;
	uint count = 0;
	for (int j = 0; j < npoints; j++) {
		float x = (float)(rng_state = rand_xorshift(rng_state)) * RAND_NORM;
		float y = (float)(rng_state = rand_xorshift(rng_state)) * RAND_NORM;
		// check if inside quarter unit circle
		if (x * x + y * y < 1.0f) {
			count++;
		}
	}
	out[i] = count;
}