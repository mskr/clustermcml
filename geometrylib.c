#include "geometrylib.h"

#define Real float
#define Real2 float2
#define Real3 float3


/*
* Distance of p to plane with origin o and normalized normal n
*/
Real calcPointPlaneDistance(Real3 p, Real3 o, Real3 n) {
    Real a = dot(o, n);
    a -= dot(n, p);
    a /= dot(n, n);
    return fabs(a);
}


/*
* An intersection routine for rays and planes
*/
Real intersectPlane(Real3 pos, Real3 dir, Real3 middle, Real3 normal) {
    Real a = dot(dir, normal);
    if (a >= 0.0) return -1.0; // facing away
    Real b = dot(middle - pos, normal);
    if (b >= 0.0) return -1.0; // behind or on plane
    return b / a;
}


/**
* Get the point on a line that minimizes the distance to another point
*/
Real3 projectPointToLine(struct Line3 line, Real3 point) {
    const Real3 s = line.start;
    const Real3 e = line.end;
    const Real num = dot(point-s, e-s);
    const Real denom = dot(e-s, e-s);
    const Real t = num / denom;
    return s + t*(e-s);
}


/**
* Get point on plane that minimizes distance to another point
* http://immersivemath.com/ila/ch03_dotproduct/ch03.html#ex_dp_ortho_proj_onto_plane
*/
Real3 projectPointToPlane(struct Plane3 plane, Real3 point) {
    Real3 v = point - plane.middle;
    Real3 proj = (dot(v, plane.normal) / pow(length(plane.normal), 2)) * plane.normal;
    return v - proj;
}


/**
* Get point rotated by angle around vector (counter-clockwise).
* Basically application of rotation matrix:
* https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
*/
Real3 rotatePointAroundVector(float angle, Real3 vector, Real3 point) {
    float cs = cos(angle), ics = 1.0 - cs, sn = sin(angle);
    Real3 u = normalize(vector);
    return (Real3)(
        point.x * (cs+u.x*u.x*ics) + point.y * (u.x*u.y*ics-u.z*sn) + point.z * (u.x*u.z*ics+u.y*sn),
        point.x * (u.y*u.x*ics+u.z*sn) + point.y * (cs+u.y*u.y*ics) + point.z * (u.y*u.z*ics-u.x*sn),
        point.x * (u.z*u.x*ics-u.y*sn) + point.y * (u.z*u.y*ics+u.x*sn) + point.z * (cs+u.z*u.z*ics));
}


/**
* Get point rotated 90 degrees around vector (counter-clockwise).
* Basically application of rotation matrix with sines and cosines canceled out.
* https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
*/
Real3 rotate90PointAroundVector(Real3 vector, Real3 point) {
    Real3 u = normalize(vector);
    return (Real3)(
        point.x * (u.x*u.x) + point.y * (u.x*u.y-u.z) + point.z * (u.x*u.z+u.y),
        point.x * (u.y*u.x+u.z) + point.y * (u.y*u.y) + point.z * (u.y*u.z-u.x),
        point.x * (u.z*u.x-u.y) + point.y * (u.z*u.y+u.x) + point.z * (u.z*u.z));
}


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
* An intersection routine for rays and cones
*
* David Eberly, Geometric Tools, Redmond WA 98052
* Copyright (c) 1998-2018
* Distributed under the Boost Software License, Version 1.0.
* http://www.boost.org/LICENSE_1_0.txt
* http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
* File Version: 3.0.2 (2018/10/05)
* Ported to OpenCL C by Marius Kircher (2018/11/20)
*
* https://www.geometrictools.com/Documentation/IntersectionLineCone.pdf
*/
Real intersectCone(Real3 lineOrigin, Real3 lineDirection, struct Cone3 cone, Real3* outNormal) {
    struct RayConeIntersectionResult result;
    // The cone has vertex V, unit-length axis direction D, angle theta in
    // (0,pi/2), and height h in (0,+infinity).  The line is P + t*U, where U
    // is a unit-length direction vector.  Define g = cos(theta).  The cone
    // is represented by
    //   (X-V)^T * (D*D^T - g^2*I) * (X-V) = 0,  0 <= dot(D,X-V) <= h
    // The first equation defines a double-sided cone.  The first inequality
    // in the second equation limits this to a single-sided cone containing
    // the ray V + s*D with s >= 0.  We will call this the 'positive cone'.
    // The single-sided cone containing ray V + s * t with s <= 0 is called
    // the 'negative cone'.  The double-sided cone is the union of the
    // positive cone and negative cone.  The second inequality in the second
    // equation limits the single-sided cone to the region bounded by the
    // height.  Setting X(t) = P + t*U, the equations are
    //   c2*t^2 + 2*c1*t + c0 = 0,  0 <= dot(D,U)*t + dot(D,P-V) <= h
    // where
    //   c2 = dot(D,U)^2 - g^2
    //   c1 = dot(D,U)*dot(D,P-V) - g^2*dot(U,P-V)
    //   c0 = dot(D,P-V)^2 - g^2*dot(P-V,P-V)
    // The following code computes the t-interval that satisfies the quadratic
    // equation subject to the linear inequality constraints.
    
    Real t;
    
    // Handle degenerate case, when there is no cone direction.
    // We interpret the cone for our use case as an annulus in the xy plane.
    // Cap and height are interpreted as inner and outer radii of the annulus.
    if (cone.ray.direction.x==0.0&&cone.ray.direction.y==0.0&&cone.ray.direction.z==0.0) {
        const Real3 normal = normalize((Real3)(0,0,-sign(lineDirection.z)));
        t = intersectPlane(lineOrigin, lineDirection, cone.ray.origin, normal);
        if (t > 0.0) {
            Real3 p = lineOrigin + t * lineDirection;
            Real r = length(p.xy - cone.ray.origin.xy);
            if (r >= cone.cap && r <= cone.height) {
                *outNormal = normal;
                return t;
            }
        }
        return -1.0;
    }

    Real3 PmV = lineOrigin - cone.ray.origin;
    Real DdU = dot(cone.ray.direction, lineDirection);
    Real DdPmV = dot(cone.ray.direction, PmV);
    Real UdPmV = dot(lineDirection, PmV);
    Real PmVdPmV = dot(PmV, PmV);
    Real cosAngle = cos(cone.angle);
    Real cosAngleSqr = cosAngle * cosAngle;
    Real c2 = DdU * DdU - cosAngleSqr;
    Real c1 = DdU * DdPmV - cosAngleSqr * UdPmV;
    Real c0 = DdPmV * DdPmV - cosAngleSqr * PmVdPmV;

    if (c2 != (Real)(0))
    {
        Real discr = c1 * c1 - c0 * c2;
        if (discr < (Real)(0))
        {
            // The quadratic has no real-valued roots.  The line does not
            // intersect the double-sided cone.
            result.type = 0;
            return -1.0;
        }
        else if (discr > (Real)(0))
        {
            // The quadratic has two distinct real-valued roots.  However, one
            // or both of them might intersect the negative cone.  We are
            // interested only in those intersections with the positive cone.
            Real root = sqrt(discr);
            Real invC2 = ((Real)(1)) / c2;
            int numParameters = 0;

            t = (-c1 - root) * invC2;
            if (DdU * t + DdPmV >= (Real)(0))
            {
                result.parameter[numParameters++] = t;
            }

            t = (-c1 + root) * invC2;
            if (DdU * t + DdPmV >= (Real)(0))
            {
                result.parameter[numParameters++] = t;
            }

            if (numParameters == 2)
            {
                // The line intersects the positive cone in two distinct
                // points.
                result.type = 2;
                if (result.parameter[0] > result.parameter[1])
                {
                    Real tmp = result.parameter[0];
                    result.parameter[0] = result.parameter[1];
                    result.parameter[1] = tmp;
                }
            }
            
            else if (numParameters == 1)
            {
                // The line intersects the positive cone in a single point and
                // the negative cone in a single point.  We report only the
                // intersection with the positive cone.
                if (DdU > (Real)(0))
                {
                    // Line enters positive cone at t==parameter[0]
                    // and keeps intersecting until t==INFINITY,
                    // in case there is no cone height constraint.
                    result.type = 3;
                    result.parameter[1] = INFINITY;
                }
                else
                {
                    // Line comes from t==-INFINITY and enters positive
                    // cone at t==parameter[0] (assuming no height constraint).
                    result.type = 4;
                    result.parameter[1] = result.parameter[0];
                    result.parameter[0] = -INFINITY;
                }
            }
            else
            {
                // The line intersects the negative cone in two distinct
                // points, but we are interested only in the intersections
                // with the positive cone.
                result.type = 0;
                return -1.0;
            }
        }
        else  // discr == 0
        {
            // One repeated real root; the line is tangent to the double-sided
            // cone at a single point.  Report only the point if it is on the
            // positive cone.
            t = -c1 / c2;
            if (DdU * t + DdPmV >= (Real)(0))
            {
                result.type = 1;
                result.parameter[0] = t;
                result.parameter[1] = t;
            }
            else
            {
                result.type = 0;
                return -1.0;
            }
        }
    }
    else if (c1 != (Real)(0))
    {
        // c2 = 0, c1 != 0; U is a direction vector on the cone boundary
        t = -((Real)(0.5))*c0 / c1;
        if (DdU * t + DdPmV >= (Real)(0))
        {
            // The line intersects the positive cone and the ray of
            // intersection is interior to the positive cone.
            if (DdU > (Real)(0))
            {
                result.type = 3;
                result.parameter[0] = t;
                result.parameter[1] = INFINITY;
            }
            else
            {
                result.type = 4;
                result.parameter[0] = -INFINITY;
                result.parameter[1] = t;
            }
        }
        else
        {
            // The line intersects the negative cone and the ray of
            // intersection is interior to the positive cone.
            result.type = 0;
            return -1.0;
        }
    }
    else if (c0 != (Real)(0))
    {
        // c2 = c1 = 0, c0 != 0.  Cross(D,U) is perpendicular to Cross(P-V,U)
        result.type = 0;
        return -1.0;
    }
    else
    {
        // c2 = c1 = c0 = 0; the line is on the cone boundary.
        result.type = 5;
        result.parameter[0] = -INFINITY;
        result.parameter[1] = +INFINITY;
    }

    // Post processing for bounded cones (cap,height)
    if (cone.height < INFINITY)
    {
        if (DdU != (Real)(0))
        {
            // Ignore intersections outside
            if (result.type > 0) {
                Real3 p0 = lineOrigin + result.parameter[0] * lineDirection;
                Real3 p1 = lineOrigin + result.parameter[1] * lineDirection;
                Real d0 = calcPointPlaneDistance(p0, cone.ray.origin, cone.ray.direction);
                Real d1 = calcPointPlaneDistance(p1, cone.ray.origin, cone.ray.direction);
                bool p0outside = (d0 <= cone.cap || d0 >= cone.height);
                bool p1outside = (d1 <= cone.cap || d1 >= cone.height);
                if (p0outside && p1outside) result.type = 0;
                else if (p0outside) result.parameter[0] = result.parameter[1];
                else if (p1outside) result.parameter[1] = result.parameter[0];
            }
        }
        else if (result.type > 0)
        {
            if (DdPmV > cone.height || DdPmV < cone.cap)
            {
                result.type = 0;
            }
        }
    }

    if (result.type > 0) {
        // Extract smallest path length to intersection
        const Real t = result.parameter[0];

        // Calc normal
        const Real3 p = lineOrigin + t*lineDirection;
        const Real3 grad = p - cone.ray.origin;
        const struct Plane3 basePlane = { cone.ray.origin, cone.ray.direction };
        const Real3 basePoint = projectPointToPlane(basePlane, p);
        const Real3 tang = rotate90PointAroundVector(cone.ray.direction, basePoint - cone.ray.origin);
        Real3 normal = cross(normalize(grad), normalize(tang));
        if (dot(normal, lineDirection) > 0.0) normal *= (Real3)(-1.0); // let normal face ray
        *outNormal = normal;

        // Return smallest path length to intersection
        return t;
    } else {
        return -1.0;
    }
}


/**
* Get i-th cone in order from inner to outer cones.
* The (capped) cone shapes arise when rotating heightfield in xy plane.
* 
* To get a cone, heighfield is sampled at i and i+1.
* Cone degenerates to annulus in xy plane, if both height samples are equal.
*
* For indices >= number of samples, annulus with infinite outer edge is returned.
* Buffer overflow is prevented by taking the last sample twice.
*/
struct Cone3 getConeAtIndex(int i, struct RHeightfield hfield) {

    // Get 2 nearest samples
    const Real h0 = hfield.heights[min(i, BOUNDARY_SAMPLES-1)];
    const Real h1 = hfield.heights[min(i+1, BOUNDARY_SAMPLES-1)];

    // Apply height
    const Real z0 = hfield.center.z - h0;
    const Real z1 = hfield.center.z - h1;

    // Calc radial offsets
    Real r0 = 0;
    for (int j = 0; j < min(i, BOUNDARY_SAMPLES-1); j++)
        r0 += hfield.spacings[j];
    Real r1 = i < BOUNDARY_SAMPLES ? r0 + hfield.spacings[i] : INFINITY;

    // Degenerate case: cone becomes annulus in xy plane
    if (h0 == h1) {
        struct Ray3 mainAxis = {
            /*origin:*/ (Real3)(hfield.center.xy, z0),
            /*normal:*/ (Real3)(0) }; // degenerate normal
        struct Cone3 annulus = { mainAxis,
             /*angle:*/  0, 
             /*height:*/ r1,
             /*cap:*/    r0 };
        return annulus;
    }

    // Linear extrapolation to get height of tip
    const Real slope = (z1-z0)/(r1-r0);
    const Real b = z0 - slope*r0;
    const Real coneTip = slope*length(hfield.center.xy) + b;

    // Construct capped cone
    struct Ray3 mainAxis = { 
        /*origin:*/ (Real3)(hfield.center.xy, coneTip), 
        /*normal:*/ h1 > h0 ? (Real3)(0,0,-1) : (Real3)(0,0,1) };
    struct Cone3 cone = { mainAxis,
         /*angle:*/  atan2(fabs(r1-r0), fabs(h1-h0)), 
         /*height:*/ max(fabs(coneTip-z0), fabs(coneTip-z1)),
         /*cap:*/    min(fabs(coneTip-z0), fabs(coneTip-z1)) };
    return cone;
}


/**
* Get index of cone from position in radial heightfield.
* The (capped) cone shapes arise when rotating heightfield in xy plane.
*
* Index starts at 0 (most inner cone).
*
* Maximum index that can be returned equals number of samples,
* i.e. is out of bounds. Should be used to handle border case.
*/
int getConeIndexFromPosition(Real2 xy, struct RHeightfield hfield) {

    // Calc p in heightmap coordinate system
    Real2 p = xy - hfield.center.xy;

    // Calc radial offset
    Real r = length(p);

    // Count spacings that fit in r
    // (clamp to number of samples)
    int i = 0;
    Real s = hfield.spacings[0];
    while (s < r) {
        i++;
        if (i < BOUNDARY_SAMPLES)
            s += hfield.spacings[i];
        else break;
    }
    return i;
}


/**
* Find intersection of line with radial heightfield (accurate)
*/
Real intersectHeightfield(struct Line3 line, struct RHeightfield hfield, Real3* outNormal) {
    // Algorithm outline:
    // iterate over all cones in the path of the line and take the closest hit
    //  - to make sure not to miss any cones, get cone under the closest point on the line to heightfield center
    //  - by assigning an ascending index from inner to outer cones, all in-between cones can easily be iterated

    const int startIndex = getConeIndexFromPosition(line.start.xy, hfield);
    const int endIndex = getConeIndexFromPosition(line.end.xy, hfield);
    const int farestIndex = max(startIndex, endIndex);
    const Real3 closestPoint = projectPointToLine(line, hfield.center);
    const int closestIndex = getConeIndexFromPosition(closestPoint.xy, hfield); 

    Real pathLenToIntersection = -1.0f;
    
    for (int i = closestIndex; i <= farestIndex; i++) {
        const struct Cone3 cone = getConeAtIndex(i, hfield);
        const Real3 lineVec = line.end - line.start;
        const Real3 rayDir = normalize(lineVec);
        Real3 normal = (Real3)(0);
        const Real t = intersectCone(line.start, rayDir, cone, &normal);
        if (t > 0.0 && t <= length(lineVec)) {
            if (t < pathLenToIntersection || pathLenToIntersection < 0.0) {
                pathLenToIntersection = t;
                *outNormal = normal;
            }
        }
    }

    return pathLenToIntersection;
}


#undef Real
#undef Real2
#undef Real3