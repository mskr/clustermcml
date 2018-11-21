#define Real float
#define Real3 float3
#define Real2 float2
#undef INFINITY
#define INFINITY 3.402823e+38


/**
*
*/
typedef struct {
    Real3 origin, direction;
} Ray3;


/**
*
*/
typedef struct {
    // The cone vertex is the ray origin and the cone axis direction is the
    // ray direction.  The direction must be unit length.  The angle must be
    // in (0,pi/2).  The height must be in (0,+infinity), where +infinity is
    // INFINITY.
    Ray3 ray;
    Real angle;
    Real height;
    Real cap; // to cap the pointy end, set this to a value in (0,height)
} Cone3;


/**
*
*/
typedef struct {
    // Because the intersection of line and cone with infinite height
    // h > 0 can be a ray or a line, we use a 'type' value that allows
    // you to decide how to interpret the parameter[] and point[] values.
    //   type  intersect  valid data
    //   0     none       none
    //   1     point      parameter[0] = parameter[1], finite
    //                    point[0] = point[1]
    //   2     segment    parameter[0] < parameter[1], finite
    //                    point[0,1] valid
    //   3     ray        parameter[0] finite, parameter[1] maxReal
    //                    point[0] = rayOrigin, point[1] = lineDirection
    //   4     ray        parameter[0] -maxReal, parameter[1] finite
    //                    point[0] = rayOrigin, point[1] = -lineDirection
    //   5     line       parameter[0] -maxReal, parameter[1] maxReal,
    //                    point[0] = lineOrigin, point[1] = lineDirection
    // If the cone height h is finite, only types 0, 1, or 2 can occur.
    int type;
    Real parameter[2];  // Relative to incoming line.
} RayConeIntersectionResult;

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
    if (a > -1e-6) return -1.0; // facing away
    Real b = dot(middle - pos, normal);
    if (b > -1e-6) return -1.0; // behind plane
    return b / a;
}

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
Real intersectCone(Real3 lineOrigin, Real3 lineDirection, Cone3 cone) {
    RayConeIntersectionResult result;
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
    // We interpret the cone for our use case as a annulus in the xy plane.
    // Cap and height are interpreted as inner and outer radii of the annulus.
    if (cone.ray.direction.x==0.0&&cone.ray.direction.y==0.0&&cone.ray.direction.z==0.0) {
        t = intersectPlane(lineOrigin, lineDirection, cone.ray.origin, (Real3)(0,0,-sign(lineDirection.z)));
        Real3 p = lineOrigin + t * lineDirection;
        Real r = length(p.xy - cone.ray.origin.xy);
        if (r >= cone.cap && r <= cone.height) {
            return t;
        } else {
            return -1.0;
        }
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

    if (cone.height < INFINITY)
    {
        if (DdU != (Real)(0))
        {
            // Ignore intersections outside (cap,height)
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
        return result.parameter[0];
    } else {
        return -1.0;
    }
}

#define SAMPLES 10
#define WIDTH 1

/**
* Return cartesian z value from radial heightfield at cartesian position xy.
* Also output normal that is oriented in -z direction.
* Heights are interpreted in -z direction since +z goes down into material.
*/
float readRadialHeightfield(Real2 pos, Real3 center, float heightfield[SAMPLES], Real3* outNormal) {
    float res = (Real)(WIDTH) / (Real)(SAMPLES);
    // Calc p in heightmap coordinate system
    Real2 p = pos.xy - center.xy;
    // Calc radial offset
    float r = length(p);
    // Calc array offset using sampling resolution
    float x = r / res;
    int i = (int)(x);
    // Get 2 nearest samples
    float h0 = heightfield[min(i, SAMPLES-1)];
    float h1 = heightfield[min(i+1, SAMPLES-1)];
    // Calc distances from the 2 samples
    float fl = 0.0; float f = fract(x, &fl); float rf = 1.0f-f;
    // Calc gradient from p1 to p0,
    // so that p0 lies on the nearest sample towards center
    // and p1 lies on the next outer sample
    Real2 p0 = (Real2)(0, 0);
    Real2 p1 = (Real2)(res, 0);
    if (x > 0.0) {
        p0 = p*(1.0f-f/x);
        p1 = p*(1.0f+rf/x);
    }
    Real3 gradient = normalize((Real3)(p0, -h0) - (Real3)(p1, -h1));
    // Calc clockwise tangent
    Real3 tangent = normalize((Real3)(p.y, -p.x, 0));
    // Calc normal from gradient and tangent
    // note: cross(a,b) forms right-handed system, where a==thumb
    // then translate normal back into cartesian coordinates
    *outNormal = normalize(cross(tangent, gradient) + center);
    // Linear interpolation
    return center.z - mix(h0, h1, f);
}

/**
* Find intersection of line with radial heightfield (inaccurate)
*/
float intersectHeightfield1(Real3 pos, Real3 dir, float len, Real3 center, float heightfield[SAMPLES], Real3* outNormal) {
    // Coarse raymarching search
    float D = min(0.02f, len);
    Real3 p = pos;
    float lastDz = p.z - readRadialHeightfield(p.xy, center, heightfield, outNormal);
    float d = D;
    Real3 ds = D * dir;
    while (d <= len) {
        p += ds;
        float dz = p.z - readRadialHeightfield(p.xy, center, heightfield, outNormal);
        // Detect sign change
        if ((lastDz > 0.0 && dz < 0.0) || (lastDz < 0.0 && dz > 0.0)) {
            // if ray came from +z, i.e. (dz < 0), normal will point down
            if (dz < 0.0) *outNormal *= -1.0f;
            // Path to intersection
            return d;
        }
        lastDz = dz;
        d += D;
    }
    return -1.f;
}

/**
* Construct the cone at given radius of the radial heightfield
*/
Cone3 getConeAt(float radius, Real3 center, float heightfield[SAMPLES]) {
    float res = (Real)(WIDTH) / (Real)(SAMPLES);
    // Calc array offset using sampling resolution
    float x = radius / res;
    int i = (int)(x);
    // Get 2 nearest samples
    float h0 = heightfield[min(i, SAMPLES-1)];
    float h1 = heightfield[min(i+1, SAMPLES-1)];
    // Calc distances from the 2 samples
    float fl = 0.0; float f = fract(x, &fl); float rf = 1.0f-f;
    // Linear interpolation
    float h = mix(h0, h1, f);
    // Apply height in -z direction
    float z = center.z - h;
    float z0 = center.z - h0;
    float z1 = center.z - h1;
    // Linearly extrapolate z to center
    float slope = (z1-z0)/res;
    float b = z - (slope*radius);
    float coneTip = slope * length(center.xy) + b;
    // Construct the cone
    return (Cone3){ (Ray3){ (Real3)(center.xy, coneTip), h1==h0 ? (Real3)(0) : h1>h0 ? (Real3)(0,0,-1):(Real3)(0,0,1) },
                 /*angle:*/  atan2(res, fabs(h1-h0)), 
                 /*height:*/ h1==h0 ? (x+rf)*res: max(fabs(coneTip-z0), fabs(coneTip-z1)),
                 /*cap:*/    h1==h0 ? (x-f)*res : min(fabs(coneTip-z0), fabs(coneTip-z1)) };
}

/**
* Find intersection of line with radial heightfield (accurate)
*/
float intersectHeightfield2(Real3 pos, Real3 dir, float len, Real3 center, float heightfield[SAMPLES], Real3* outNormal) {
    // iterate over all cones in the path of the ray and take the first hit
    float res = (Real)(WIDTH) / (Real)(SAMPLES);
    float s = 0.0;
    while (s < len) {
        Real3 p = pos + s * dir;
        Cone3 cone = getConeAt(length(p.xy - center.xy), center, heightfield);
        Real t = intersectCone(pos, dir, cone);
        if (t > 0) {
            
            //TODO get normal
            
            return t;
        }
        s += res;
    }
    return -1.f;
}