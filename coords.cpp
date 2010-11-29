#include "common.h"
#include "coords.h"


Vector3 createSomePerpendicular(const Vector3& v)
{
	Vector3 perp;
	Vector3 par = Vector3(v).normalize();

	Float eps = v.norm() * kMachineEpsilon;

	if (fabs(v.y()) > eps || fabs(v.z()) > eps)
		perp = Vector3(1., 0., 0.);
	else
		perp = Vector3(0., 1., 0.);
	
	perp -= par*(perp*par);

	return perp;
}

Vector3 createSomeDeviantVector(const Vector3& v, const Angle& a)
{
	Vector3 par = Vector3(v).normalize();
	Vector3 perp = createSomePerpendicular(v);

	return par*a.costheta + perp*a.sintheta;
}

Matrix3 createTransformMatrix(const Vector3& v1, const Vector3& v2, const Vector3& v3,
							  const Vector3& u1, const Vector3& u2, const Vector3& u3)
{
	Matrix3 m;

	m(0,0) = v1*u1;
	m(0,1) = v1*u2;
	m(0,2) = v1*u3;

	m(1,0) = v2*u1;
	m(1,1) = v2*u2;
	m(1,2) = v2*u3;

	m(2,0) = v3*u1;
	m(2,1) = v3*u2;
	m(2,2) = v3*u3;

	return m;
}

Matrix3 createTransformMatrix(const Vector3& u1, const Vector3& u2, const Vector3& u3)
{
	Matrix3 m;

	m(0,0) = u1.x();
	m(1,0) = u2.x();
	m(2,0) = u3.x();

	m(0,1) = u1.y();
	m(1,1) = u2.y();
	m(2,1) = u3.y();

	m(0,2) = u1.z();
	m(1,2) = u2.z();
	m(2,2) = u3.z();

	return m;
}



