// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.2 (2018/10/05)
// Adapted to fit OpenCL C by Marius Kircher (2018/11/08)

#define Real float
#define Real3 float3
#define INFINITY 3.402823e+38


/**
*
*/
struct Ray {
    Real3 origin, direction;
};


/**
*
*/
struct Cone3 {
    // The cone vertex is the ray origin and the cone axis direction is the
    // ray direction.  The direction must be unit length.  The angle must be
    // in (0,pi/2).  The height must be in (0,+infinity), where +infinity is
    // INFINITY.
    Ray ray;
    Real angle;
    Real height;
};


/**
*
*/
struct Result {
    bool intersect;

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
    Real3 point[2];
};


/**
*
*/
struct IntervalIntervalResult {
    bool intersect;

    // Static queries (no motion of intervals over time).  The number of
    // number of intersections is 0 (no overlap), 1 (intervals are just
    // touching), or 2 (intervals overlap in an interval).  If 'intersect'
    // is false, numIntersections is 0 and 'overlap' is set to
    // [maxReal,-maxReal].  If 'intersect' is true, numIntersections is
    // 1 or 2.  When 1, 'overlap' is set to [x,x], which is degenerate and
    // represents the single intersection point x.  When 2, 'overlap' is
    // the interval of intersection.
    int numIntersections;
    Real overlap[2];

    // Dynamic queries (intervals moving with constant speeds).  If
    // 'intersect' is true, the contact times are valid and
    //     0 <= firstTime <= lastTime,  firstTime <= maxTime
    // If 'intersect' is false, there are two cases reported.  If the
    // intervals will intersect at firstTime > maxTime, the contact times
    // are reported just as when 'intersect' is true.  However, if the
    // intervals will not intersect, then firstTime = maxReal and
    // lastTime = -maxReal.
    Real firstTime, lastTime;
};


/**
* Intersection query for two intervals (1-dimensional query).
*
* The intervals are [u0,u1] and [v0,v1], where u0 <= u1 and v0 <= v1, and
* where the endpoints are any finite floating-point numbers.  Degenerate
* intervals are allowed (u0 = u1 or v0 = v1).  The query does not perform
* validation on the input intervals.
*
* https://www.geometrictools.com/Documentation/IntersectionLine2Circle2.pdf
*/
IntervalIntervalResult findIntersectionIntervalInterval(Real* interval0, Real* interval1) {
    IntervalIntervalResult result;
    result.firstTime = INFINITY;
    result.lastTime = -INFINITY;

    if (interval0[1] < interval1[0] || interval0[0] > interval1[1])
    {
        result.numIntersections = 0;
        result.overlap[0] = INFINITY;
        result.overlap[1] = -INFINITY;
    }
    else if (interval0[1] > interval1[0])
    {
        if (interval0[0] < interval1[1])
        {
            result.numIntersections = 2;
            result.overlap[0] =
                (interval0[0] < interval1[0] ? interval1[0] : interval0[0]);
            result.overlap[1] =
                (interval0[1] > interval1[1] ? interval1[1] : interval0[1]);
            if (result.overlap[0] == result.overlap[1])
            {
                result.numIntersections = 1;
            }
        }
        else  // interval0[0] == interval1[1]
        {
            result.numIntersections = 1;
            result.overlap[0] = interval0[0];
            result.overlap[1] = result.overlap[0];
        }
    }
    else  // interval0[1] == interval1[0]
    {
        result.numIntersections = 1;
        result.overlap[0] = interval0[1];
        result.overlap[1] = result.overlap[0];
    }

    result.intersect = (result.numIntersections > 0);
    return result;
}


/**
* https://www.geometrictools.com/Documentation/IntersectionLineCone.pdf
*/
Result intersectCone(Real3 lineOrigin, Real3 lineDirection, Cone3 cone) {
    Result result;
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
    Real t;

    if (c2 != (Real)0)
    {
        Real discr = c1 * c1 - c0 * c2;
        if (discr < (Real)0)
        {
            // The quadratic has no real-valued roots.  The line does not
            // intersect the double-sided cone.
            result.intersect = false;
            result.type = 0;
            return result;
        }
        else if (discr > (Real)0)
        {
            // The quadratic has two distinct real-valued roots.  However, one
            // or both of them might intersect the negative cone.  We are
            // interested only in those intersections with the positive cone.
            Real root = sqrt(discr);
            Real invC2 = ((Real)1) / c2;
            int numParameters = 0;

            t = (-c1 - root) * invC2;
            if (DdU * t + DdPmV >= (Real)0)
            {
                result.parameter[numParameters++] = t;
            }

            t = (-c1 + root) * invC2;
            if (DdU * t + DdPmV >= (Real)0)
            {
                result.parameter[numParameters++] = t;
            }

            if (numParameters == 2)
            {
                // The line intersects the positive cone in two distinct
                // points.
                result.intersect = true;
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
                result.intersect = true;
                if (DdU > (Real)0)
                {
                    result.type = 3;
                    result.parameter[1] = INFINITY;
                }
                else
                {
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
                result.intersect = false;
                result.type = 0;
                return result;
            }
        }
        else  // discr == 0
        {
            // One repeated real root; the line is tangent to the double-sided
            // cone at a single point.  Report only the point if it is on the
            // positive cone.
            t = -c1 / c2;
            if (DdU * t + DdPmV >= (Real)0)
            {
                result.intersect = true;
                result.type = 1;
                result.parameter[0] = t;
                result.parameter[1] = t;
            }
            else
            {
                result.intersect = false;
                result.type = 0;
                return result;
            }
        }
    }
    else if (c1 != (Real)0)
    {
        // c2 = 0, c1 != 0; U is a direction vector on the cone boundary
        t = -((Real)0.5)*c0 / c1;
        if (DdU * t + DdPmV >= (Real)0)
        {
            // The line intersects the positive cone and the ray of
            // intersection is interior to the positive cone.
            result.intersect = true;
            if (DdU > (Real)0)
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
            result.intersect = false;
            result.type = 0;
            return result;
        }
    }
    else if (c0 != (Real)0)
    {
        // c2 = c1 = 0, c0 != 0.  Cross(D,U) is perpendicular to Cross(P-V,U)
        result.intersect = false;
        result.type = 0;
        return result;
    }
    else
    {
        // c2 = c1 = c0 = 0; the line is on the cone boundary.
        result.intersect = true;
        result.type = 5;
        result.parameter[0] = -INFINITY;
        result.parameter[1] = +INFINITY;
    }

    if (cone.height < INFINITY)
    {
        if (DdU != (Real)0)
        {
            // Clamp the intersection to the height of the cone.
            Real invDdU = ((Real)1) / DdU;
            Real[2] hInterval;
            if (DdU >(Real)0)
            {
                hInterval[0] = -DdPmV * invDdU;
                hInterval[1] = (cone.height - DdPmV) * invDdU;
            }
            else // (DdU < (Real)0)
            {
                hInterval[0] = (cone.height - DdPmV) * invDdU;
                hInterval[1] = -DdPmV * invDdU;
            }

            IntervalIntervalResult iiResult = findIntersectionIntervalInterval(result.parameter, hInterval);
            result.intersect = (iiResult.numIntersections > 0);
            result.type = iiResult.numIntersections;
            result.parameter = iiResult.overlap;
        }
        else if (result.intersect)
        {
            if (DdPmV > cone.height)
            {
                result.intersect = false;
                result.type = 0;
            }
        }
    }

    return result;
}
