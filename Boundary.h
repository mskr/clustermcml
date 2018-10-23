#define BOUNDARY_SAMPLES 100
#define BOUNDARY_WIDTH 5 // cm

ALIGN_4BYTE(
struct Boundary {
	float z; // origin (0,0,z)
	float heightfield[BOUNDARY_SAMPLES];
});