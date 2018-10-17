#define BOUNDARY_SAMPLES 100

struct Boundary {
	float ox, oy, oz; // origin
	float heightfield[BOUNDARY_SAMPLES];
};