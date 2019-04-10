#define BOUNDARY_SAMPLES 50

#include "geometrylib.c.cpp"
#define Real float
#define Real3 vec3
#include <stdio.h>

/**
* Return index difference from current to adjacent layer in case
* there is an intersection with the boundary in the photon step.
* Also outputs collision normal and path length to intersection.
*/
int detectBoundaryCollision(int currentLayer, Line3 line, RHeightfield* boundaries, float3* normal, float* pathLenToIntersection) {
	// Find intersection with top boundary
	*pathLenToIntersection = intersectHeightfield(line, boundaries[currentLayer], normal);
	if (*pathLenToIntersection >= 0) {
		return -1;
	}
	// Find intersection with bottom boundary
	*pathLenToIntersection = intersectHeightfield(line, boundaries[currentLayer+1], normal);
	if (*pathLenToIntersection >= 0) {
		return 1;
	}
	return 0;
}


RHeightfield boundaries[1] = {
	{
		{0.0f,0.0f,0.0f},
		{0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f},
		{0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f,0.1f}
	}
};

int main() {
	int currentLayer = 0;
	Line3 line = { {-0.0127144791,-0.00177060463,9.85572115e-07}, {-0.0155298375,-0.00247456646,-0.00107741705} };
	Real3 normal;
	Real pathLenToIntersection;
	int layerChange = detectBoundaryCollision(currentLayer, line, boundaries, &normal, &pathLenToIntersection);
	printf("layerChange=%d normal=(%f,%f,%f) pathLenToIntersection=%f\n", layerChange, normal.x, normal.y, normal.z, pathLenToIntersection);
	if (layerChange == -1) {
		printf("TEST PASSED\n");
		return 0;
	} else {
		printf("TEST FAILED\n");
		return 1;
	}

	//TODO Test fails => debug this
	// same as finishedPhotonCount = 17593, iteration = 1
}