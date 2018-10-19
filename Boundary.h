#define BOUNDARY_SAMPLES 100

ALIGN_4BYTE(
struct Boundary {
	float z; // origin (0,0,z)
	float heightfield[BOUNDARY_SAMPLES];
});