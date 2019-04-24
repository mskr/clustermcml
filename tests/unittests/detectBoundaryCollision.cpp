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







float heights[10] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
float spacings[10] = {0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f,0.01f};

Boundary boundaries[2] = {
	// {
	// 	1, 0, 0, 0, RHeightfield{ // This boundary is a heightfield
	// 		{0.0f, 0.0f, 0.0f}, // origin at 0
	// 		0, 5, // first half of the sample arrays
	// 		0, 5
	// 	}
	// },
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



int main() {

	//TODO no reflectance detected when using implicit boundary => bug in detectBoundaryCollision?


	int currentLayer = 0;
	Line3 line = { {0.0f, 0.0f, 0.05f}, {0.0f, 0.0f, -0.05f} };
	float3 normal;
	float pathLenToIntersection;

	int layerChange = detectBoundaryCollision(currentLayer, line, boundaries, heights, spacings, &normal, &pathLenToIntersection);

	printf("layerChange=%d normal=(%f,%f,%f) pathLenToIntersection=%f\n", layerChange, normal.x, normal.y, normal.z, pathLenToIntersection);

	if (layerChange == -1) {
		printf("TEST PASSED\n");
		return 0;
	} else {
		printf("TEST FAILED\n");
		return 1;
	}
}