#ifndef COORDS_CQTZQF5F

#define COORDS_CQTZQF5F

#include "vector3.h"
#include "matrix3.h"
#include "angle.h"

Vector3 createSomePerpendicular(const Vector3& v);
Vector3 createSomeDeviantVector(const Vector3& v, const Angle& a); //lol, funny function name


///matrix to transform vector from (v1,v2,v3) to (u1,u2,u3) coordinate system
Matrix3 createTransformMatrix(  const Vector3& v1, const Vector3& v2, const Vector3& v3,
								const Vector3& u1, const Vector3& u2, const Vector3& u3);

#endif /* end of include guard: COORDS_CQTZQF5F */
