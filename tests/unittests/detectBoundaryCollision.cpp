#include <stdio.h>
#include <stdint.h>

#include "Layer.h" // struct Layer
// typedef struct { float x,y,z; } float3; // remove when transpiling
#include "Boundary.h" // struct Boundary

#define __global
#define BoundaryArray __global struct Boundary*

/**
* Return index difference from current to adjacent layer in case
* there is an intersection with the boundary in the photon step.
* Also outputs collision normal and path length to intersection.
*/
// int detectBoundaryCollision_old(int currentLayer, Line3 line, RHeightfield* boundaries, float3* normal, float* pathLenToIntersection) {
// 	// Find intersection with top boundary
// 	*pathLenToIntersection = intersectHeightfield(line, boundaries[currentLayer], normal);
// 	if (*pathLenToIntersection >= 0) {
// 		return -1;
// 	}
// 	// Find intersection with bottom boundary
// 	*pathLenToIntersection = intersectHeightfield(line, boundaries[currentLayer+1], normal);
// 	if (*pathLenToIntersection >= 0) {
// 		return 1;
// 	}
// 	return 0;
// }

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







struct Result {
	float3 normal;
	float pathLenToIntersection;
	int layerChange;
};


Result test(int currentLayer, Line3 line, Boundary* boundaries, float* heights, float* spacings) {

	printf("line={ start=(%f,%f,%f) end=(%f,%f,%f) }\n", line.start[0], line.start[1], line.start[2], line.end[0], line.end[1], line.end[2]);
	printf("topBoundary=   { isHeightfield=%d", boundaries[currentLayer].isHeightfield);
	if (!boundaries[currentLayer].isHeightfield) printf(", z=%f", boundaries[currentLayer].z);
	else {
		printf(", origin=(%f,%f,%f)", boundaries[currentLayer].heightfield.center[0], boundaries[currentLayer].heightfield.center[1], boundaries[currentLayer].heightfield.center[2]);
		printf(", heights="); for(unsigned int i = boundaries[currentLayer].heightfield.i_heights; i < boundaries[currentLayer].heightfield.i_heights+boundaries[currentLayer].heightfield.n_heights; i++) printf("%f ", heights[i]);
	}
	printf("\n");
	printf("bottomBoundary={ isHeightfield=%d", boundaries[currentLayer+1].isHeightfield);
	if (!boundaries[currentLayer+1].isHeightfield) printf(", z=%f }", boundaries[currentLayer+1].z);
	else {
		printf(", origin=(%f,%f,%f)", boundaries[currentLayer+1].heightfield.center[0], boundaries[currentLayer+1].heightfield.center[1], boundaries[currentLayer+1].heightfield.center[2]);
		printf(", heights="); for(unsigned int i = boundaries[currentLayer+1].heightfield.i_heights; i < boundaries[currentLayer+1].heightfield.i_heights+boundaries[currentLayer+1].heightfield.n_heights; i++) printf("%f ", heights[i]);
	}
	printf("\n");

	Result result;
	result.layerChange = detectBoundaryCollision(currentLayer, line, boundaries, heights, spacings, &result.normal, &result.pathLenToIntersection);

	printf("layerChange=%d normal=(%f,%f,%f) pathLenToIntersection=%f\n", result.layerChange, result.normal.x, result.normal.y, result.normal.z, result.pathLenToIntersection);

	return result;
}



int main() {

	int currentLayer = 0;
	Line3 line = { {0.0f, 0.0f, 0.05f}, {100.0f, 100.0f, -0.05f} };

	float heights[10] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
	float spacings[10] = {0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f};

	Boundary boundaries1[2] = {
		{
			1, 0, 0, 0, RHeightfield{ // This boundary is a heightfield
				{0.0f, 0.0f, 0.0f}, // origin at 0
				0, 5, // first half of the sample arrays
				0, 5
			}
		},
		{
			1, 0, 0, 0, RHeightfield{ // This boundary is a heightfield
				{0.0f, 0.0f, 0.1f}, // origin at z=0.1
				5, 5, // second half of the sample arrays
				5, 5
			}
		}
	};

	Result result1 = test(currentLayer, line, boundaries1, heights, spacings);


	if (result1.layerChange == -1) {
		printf("TEST 1 PASSED\n");
	} else {
		printf("TEST 1 FAILED\n");
		return 1;
	}


	Boundary boundaries2[2] = {
		{
			0, 0.0f, 0, 0, RHeightfield{ // This boundary is implicit
				{0.0f, 0.0f, 0.0f},
				0, 5,
				0, 5
			}
		},
		{
			1, 0, 0, 0, RHeightfield{ // This boundary is a heightfield
				{0.0f, 0.0f, 0.1f}, // origin at z=0.1
				5, 5, // second half of the sample arrays
				5, 5
			}
		}
	};

	Result result2 = test(currentLayer, line, boundaries2, heights, spacings);


	if (result2.layerChange == -1) {
		printf("TEST 2 PASSED\n");
	} else {
		printf("TEST 2 FAILED\n");
		return 1;
	}

	return 0;
}