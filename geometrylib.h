#define Real float
#define Real3 float3

/**
* Ray in 3-space
*/
struct Ray3 {
    Real3 origin, direction;
};


/**
* Cone in 3-space.
* The cone has vertex V, unit-length axis direction D, angle theta in (0,pi/2).
* Also a height and cap position are defined along the axial ray.
* The cone vertex is the ray origin and the cone axis direction is the
* ray direction. The direction must be unit length. The angle must be
* in (0,pi/2). The height must be in (0,+infinity), where +infinity is INFINITY.
*/
struct Cone3 {
    struct Ray3 ray;
    Real angle;
    Real height;
    Real cap; // to cap the pointy end, set this to a value in (0,height)
};


/**
* Line in 3-space
*/
struct Line3 {
    Real3 start;
    Real3 end;
};


/**
*
*/
struct Plane3 {
    Real3 middle;
    Real3 normal;
};


/**
* Radial heightfield
*/
struct RHeightfield {
    Real3 center;
    Real heights[BOUNDARY_SAMPLES];
    Real spacings[BOUNDARY_SAMPLES];
};


/**
* Because the intersection of line and cone with infinite height
* h > 0 can be a ray or a line, we use a 'type' value that allows
* you to decide how to interpret the parameter[] and point[] values.
*   type  intersect  valid data
*   0     none       none
*   1     point      parameter[0] = parameter[1], finite
*                    point[0] = point[1]
*   2     segment    parameter[0] < parameter[1], finite
*                    point[0,1] valid
*   3     ray        parameter[0] finite, parameter[1] maxReal
*                    point[0] = rayOrigin, point[1] = lineDirection
*   4     ray        parameter[0] -maxReal, parameter[1] finite
*                    point[0] = rayOrigin, point[1] = -lineDirection
*   5     line       parameter[0] -maxReal, parameter[1] maxReal,
*                    point[0] = lineOrigin, point[1] = lineDirection
* If the cone height h is finite, only types 0, 1, or 2 can occur.
*/
struct RayConeIntersectionResult {
    int type;
    Real parameter[2];  // Relative to incoming line.
};


/**
* Find intersection of line with radial heightfield (accurate)
*/
Real intersectHeightfield(struct Line3 line, struct RHeightfield hfield, Real3* outNormal);

#undef Real
#undef Real3