#ifndef GEOMETRYLIB_H
#define GEOMETRYLIB_H

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
* Holds center, but not actual sample values,
* because they are dynamic. Instead it points
* to ranges in arrays of heights and spacings.
*/
struct RHeightfield {
    Real3 center;

    uint i_heights;
    uint n_heights;

    uint i_spacings;
    uint n_spacings;
};


/**
* Distance of p to plane with origin o and normalized normal n
*/
Real calcPointPlaneDistance(Real3 p, Real3 o, Real3 n) ;


/**
* An intersection routine for rays and planes
*/
Real intersectPlane(Real3 pos, Real3 dir, Real3 middle, Real3 normal);


/**
* Plane-line intersection implemented on top of ray-plane intersection
*/
Real intersectPlaneWithLine(struct Plane3 plane, struct Line3 line);


/**
* Get the point on a line that minimizes the distance to another point
*/
Real3 projectPointToLine(struct Line3 line, Real3 point);


/**
* Get point on plane that minimizes distance to another point
* http://immersivemath.com/ila/ch03_dotproduct/ch03.html#ex_dp_ortho_proj_onto_plane
*/
Real3 projectPointToPlane(struct Plane3 plane, Real3 point);


/**
* Get point rotated by angle around vector (counter-clockwise).
* Basically application of rotation matrix:
* https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
*/
Real3 rotatePointAroundVector(float angle, Real3 vector, Real3 point);


/**
* Get point rotated 90 degrees around vector (counter-clockwise).
* Basically application of rotation matrix with sines and cosines canceled out.
* https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
*/
Real3 rotate90PointAroundVector(Real3 vector, Real3 point);


/**
* Find intersection of line with radial heightfield (accurate)
*/
Real intersectHeightfield(struct Line3 line, struct RHeightfield hfield, __global Real* heights, __global Real* spacings, Real3* outNormal);

#undef Real
#undef Real3

#endif